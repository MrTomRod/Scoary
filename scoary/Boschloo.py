import os
import sqlite3
from typing import Optional

import pandas as pd
import numpy as np
from multiprocessing import Pool

from .simple_boschloo import simple_boschloo

root = os.path.dirname(__file__)


# root = '/home/thomas/PycharmProjects/pfh-metabolome/orthogene_analysis_resolution_tax/boschloo'


def boschloo_swap(c1r1: int, c2r1: int, c1r2: int, c2r2: int) -> str:
    """
    Two contingency tables always give the same result:

    boschloo_exact([[c1r1, c2r1], [c1r2, c2r2]])
    and
    boschloo_exact([[c2r1, c1r1], [c2r2, c1r2]])
    are equivalent.

    Compute and save only one version.
    """
    if c1r1 < c2r1:
        return f'{c1r1},{c2r1},{c1r2},{c2r2}'
    else:
        return f'{c2r1},{c1r1},{c2r2},{c1r2}'


def boschloo_swap_series(data: pd.Series):
    return boschloo_swap(*(data[x] for x in ('c1r1', 'c2r1', 'c1r2', 'c2r2')))


list_to_string = lambda l: ', '.join(f"'{e}'" for e in l)
list_to_string_bracket = lambda l: ', '.join(f"('{e}')" for e in l)


class Boschloo:
    def __init__(self, db_path: str = None, n_cpus: int = None):
        if db_path is None:
            if 'BOSCHLOO_DB' in os.environ:
                db_path = os.environ['BOSCHLOO_DB']
            else:
                db_path = os.path.expanduser('~/.cache/boschloo.db')

        self.n_cpus = n_cpus if n_cpus else os.cpu_count()
        self._db_path = db_path
        db_exists = os.path.isfile(db_path)
        self.con, self.cur = self.get_cur()
        if not db_exists:
            self.__create_db()

    def get_cur(self):
        con = sqlite3.connect(self._db_path)
        cur = con.cursor()
        return con, cur

    def __del__(self):
        self.cur.close()
        self.con.close()

    def __create_db(self):
        self.cur.execute(f'''
        CREATE TABLE boschloo (
            test text,
            pval real,
            stat real,
            PRIMARY KEY (test)
            )''')

    def get_or_create(self, c1r1: int, c2r1: int, c1r2: int, c2r2: int) -> (float, float):
        test_string = boschloo_swap(c1r1, c2r1, c1r2, c2r2)

        res = self.cur.execute(
            'SELECT pval, stat FROM boschloo WHERE test = ?',
            (test_string,)
        ).fetchone()

        if res is None:
            res = self._create(test_string)

        return res

    def _create(self, test_string: str) -> (float, float):
        pval, stat = simple_boschloo(test_string)
        self.cur.execute(f'INSERT OR IGNORE INTO boschloo VALUES (?, ?, ?)', (test_string, pval, stat))
        self.con.commit()
        return pval, stat

    def _create_many(self, test_strings: [str]):
        with Pool(self.n_cpus) as p:
            res = p.map(simple_boschloo, test_strings)

        res = [(test_string, pval, stat) for test_string, (pval, stat) in zip(test_strings, res)]

        self.cur.executemany(
            'INSERT OR IGNORE INTO boschloo VALUES (?, ?, ?)',
            res
        )
        self.con.commit()

        return pd.DataFrame(res, columns=['test', 'pval', 'stat'])

    def get_or_create_many(self, boschloo_df: pd.DataFrame, create_only: bool = False) -> Optional[pd.DataFrame]:
        assert list(boschloo_df.columns) == ['c1r1', 'c2r1', 'c1r2', 'c2r2'], f'boschloo_columns do not match! {boschloo_df.columns=}'
        boschloo_df['test'] = boschloo_df.apply(boschloo_swap_series, axis=1)
        unique_tests = boschloo_df['test'].unique()

        missing = self.cur.execute(
            f'''
            WITH mytable (test) AS ( VALUES {list_to_string_bracket(unique_tests)} )
            SELECT test
            FROM mytable
            WHERE NOT EXISTS( SELECT test FROM boschloo WHERE boschloo.test = mytable.test )
            '''

        ).fetchall()

        missing = [m[0] for m in missing]

        self._create_many(missing)

        print(f'Calculated {len(missing) / len(unique_tests):.1%} ({len(missing)} out of {len(unique_tests)}) Boschloo\'s tests')

        if create_only:
            return

        res = self.cur.execute(
            f'SELECT test, pval, stat FROM boschloo WHERE test IN ({list_to_string(unique_tests)})'
        ).fetchall()

        res = pd.DataFrame(data=res, columns=['test', 'pval', 'stat'])

        res = pd.merge(boschloo_df, res, on=['test'], how='left')

        assert not any(res.pval.isna()), f'Programming error: Failed to computes some Boschloo tests!'

        res.drop('test', axis=1, inplace=True)

        return res


def test_many():
    boschloo = Boschloo()
    boschloo_df = pd.DataFrame(
        [(np.random.randint(100), np.random.randint(100), np.random.randint(100), np.random.randint(100)) for x in range(10)],
        columns=['c1r1', 'c2r1', 'c1r2', 'c2r2']
    )
    return boschloo.get_or_create_many(boschloo_df)


# res = test_many()


def test_swap():
    import time
    boschloo = Boschloo()
    for i in range(10):
        c1r1, c2r1, c1r2, c2r2 = np.random.randint(100), np.random.randint(100), np.random.randint(100), np.random.randint(100)
        boschloo.get_or_create(c1r1, c2r1, c1r2, c2r2)
        start = time.monotonic()
        boschloo.get_or_create(c2r1, c1r1, c2r2, c1r2)
        assert (time.monotonic() - start) < 0.001, f'Failed to cache somehow'

# test_swap()
