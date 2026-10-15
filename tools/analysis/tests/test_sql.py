"""Queries name relations only for as long as they run (`pt_analysis.sql.Scope`): deriving, querying and deriving
again leaves no name on the package's connection for a later query to collide with."""

import shutil

import duckdb
import pytest

from pt_analysis import edge_stats, flows, network, pair, strata
from pt_analysis.derive import analysis_dir, derive_analysis, derived_dir
from pt_analysis.sql import Scope, leftover_names
from pt_analysis.store import CONNECTION


def query_everything(analysis):
    strata.stratified_functionals(analysis)
    flows.reassign(analysis, max_lag=100)
    edge_stats.edge_table(analysis)
    network.network(analysis)
    return pair.inspect_pair(analysis, "1", "3")


def test_derive_query_and_derive_again_leave_no_names_behind(golden_analysis):
    assert derive_analysis(golden_analysis) is True
    first = query_everything(golden_analysis)
    assert leftover_names(CONNECTION) == []
    for store in golden_analysis.runs.values():
        shutil.rmtree(derived_dir(store))
    shutil.rmtree(analysis_dir(golden_analysis))
    assert derive_analysis(golden_analysis) is True
    assert query_everything(golden_analysis) == first
    assert leftover_names(CONNECTION) == []


def test_a_scope_unregisters_its_names_even_when_its_query_fails():
    connection = duckdb.connect()
    with pytest.raises(duckdb.Error):
        with Scope(connection) as s:
            s.rows(f"SELECT missing FROM {s.name(connection.sql('SELECT 1 AS x'))}")
    assert leftover_names(connection) == []


def test_the_package_connection_does_not_read_python_locals_by_name():
    local = CONNECTION.sql("SELECT 1 AS x")  # noqa: F841  (what a replacement scan would find)
    with pytest.raises(duckdb.CatalogException):
        CONNECTION.sql("SELECT * FROM local").fetchall()
