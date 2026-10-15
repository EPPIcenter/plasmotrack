import { afterEach, describe, expect, it } from "vitest";
import { clearSelection, getSelection, select, selectOnly, subscribeSelection } from "./selection";

afterEach(() => clearSelection());

describe("the shared Infection selection", () => {
  it("starts empty", () => {
    expect([...getSelection()]).toEqual([]);
  });

  it("a click selects one Infection, replacing the selection", () => {
    select("3");
    select("1");
    expect([...getSelection()]).toEqual(["1"]);
  });

  it("a shift-click adds an Infection, and removes it when already selected", () => {
    select("3");
    select("1", { add: true });
    expect([...getSelection()].sort()).toEqual(["1", "3"]);
    select("3", { add: true });
    expect([...getSelection()]).toEqual(["1"]);
  });

  it("tells subscribers of every change until they unsubscribe", () => {
    const seen: string[][] = [];
    const unsubscribe = subscribeSelection((ids) => seen.push([...ids].sort()));
    select("2");
    select("1", { add: true });
    clearSelection();
    unsubscribe();
    select("3");
    expect(seen).toEqual([["2"], ["1", "2"], []]);
  });

  it("selects a set of Infections at once, replacing the selection (an edge's two ends)", () => {
    select("5");
    selectOnly(["1", "3"]);
    expect([...getSelection()].sort()).toEqual(["1", "3"]);
  });
});
