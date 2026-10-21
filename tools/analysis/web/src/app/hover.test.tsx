import { act, fireEvent, render, renderHook, screen } from "@testing-library/react";
import { describe, expect, it } from "vitest";
import { DataTable, type Column } from "../components/DataTable";
import { hover, isHovered, useHovered } from "./hover";

describe("the hover store", () => {
  it("holds the entity under the pointer for every chart to read, null when none", () => {
    const { result } = renderHook(() => useHovered());
    expect(result.current).toBeNull();
    act(() => hover.set({ kind: "infection", id: "pf0147" }));
    expect(result.current).toEqual({ kind: "infection", id: "pf0147" });
    expect(isHovered(result.current, { kind: "infection", id: "pf0147" })).toBe(true);
    expect(isHovered(result.current, { kind: "infection", id: "pf0001" })).toBe(false);
    act(() => hover.set(null));
    expect(result.current).toBeNull();
  });

  it("is set by a DataTable row's entity on hover and cleared on leaving the table", () => {
    interface Row { id: string }
    const columns: Column<Row>[] = [{ id: "id", header: "Infection", value: (r) => r.id }];
    render(<DataTable rows={[{ id: "a" }, { id: "b" }]} columns={columns} getRowId={(r) => r.id} hoverEntity={(r) => ({ kind: "infection", id: r.id })} />);
    fireEvent.mouseEnter(screen.getAllByRole("row")[1]!);
    expect(hover.get()).toEqual({ kind: "infection", id: "b" });
    fireEvent.mouseLeave(screen.getByRole("table"));
    expect(hover.get()).toBeNull();
  });
});
