"""Epi data: a CSV keyed by Infection id, read with DuckDB, typed by inspection and checked against the Infections
of the Analysis's Runs."""

import pytest
from conftest import open_golden

from pt_analysis.epi import EpiError, load_epi

CSV = """infection,district,observed,age
1,North,2020-01-05,4
3,North,2020-03-10,31
2,South,2020-02-01,12.5
9,South,2020-02-02,7
"""


def open_analysis(tmp_path, csv=CSV, name="epi.csv", **sections):
    if csv is not None:
        (tmp_path / name).write_text(csv)
    return open_golden(tmp_path, epi=name, **sections)


def test_columns_are_typed_by_inspection(tmp_path):
    epi = load_epi(open_analysis(tmp_path))
    assert epi.key == "infection"
    assert epi.columns == {"district": "categorical", "observed": "date", "age": "numeric"}


def test_ids_without_an_infection_are_reported_and_infections_without_a_row_are_listed(tmp_path):
    epi = load_epi(open_analysis(tmp_path, csv=CSV.replace("2,South,2020-02-01,12.5\n", "")))
    assert epi.unmatched == ["9"]
    assert epi.missing == ["2"]


def test_each_infection_has_its_value_of_a_column(tmp_path):
    epi = load_epi(open_analysis(tmp_path))
    assert epi.values("district") == {"1": "North", "2": "South", "3": "North"}
    assert epi.values("observed") == {"1": "2020-01-05", "2": "2020-02-01", "3": "2020-03-10"}
    assert epi.values("age") == {"1": "4", "2": "12.5", "3": "31"}


def test_without_an_infection_column_the_first_column_is_the_key(tmp_path):
    epi = load_epi(open_analysis(tmp_path, csv="id,travel\n1,yes\n2,no\n3,\n"))
    assert epi.key == "id"
    assert epi.columns == {"travel": "categorical"}
    assert epi.values("travel") == {"1": "yes", "2": "no", "3": None}


def test_ids_are_read_as_text(tmp_path):
    epi = load_epi(open_analysis(tmp_path, csv="infection,x\n01,a\n1,b\n2,c\n3,d\n"))
    assert epi.unmatched == ["01"]


def test_a_missing_column_is_a_named_error(tmp_path):
    epi = load_epi(open_analysis(tmp_path))
    with pytest.raises(EpiError, match="no column 'travel'.*district, observed, age"):
        epi.values("travel")


def test_an_analysis_without_epi_data_has_none(golden_analysis):
    assert load_epi(golden_analysis) is None


def test_a_missing_file_or_a_repeated_id_is_an_error(tmp_path):
    with pytest.raises(EpiError, match="not found"):
        load_epi(open_analysis(tmp_path, csv=None))
    with pytest.raises(EpiError, match="more than one row.*'1'"):
        load_epi(open_analysis(tmp_path, csv="infection,x\n1,a\n1,b\n"))


def test_stratum_values_come_in_the_columns_own_order_with_no_value_last(tmp_path):
    epi = load_epi(open_analysis(tmp_path, csv="infection,age,when,group\n1,10,2020-02-01,b\n2,9,2020-01-15,\n3,,2019-12-31,a\n"))
    assert epi.stratum_values("age") == ["9", "10", None]
    assert epi.stratum_values("when") == ["2019-12-31", "2020-01-15", "2020-02-01"]
    assert epi.stratum_values("group") == ["a", "b", None]
