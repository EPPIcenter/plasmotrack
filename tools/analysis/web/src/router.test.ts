import { describe, expect, it } from "vitest";
import { formatRoute, parseRoute } from "./router";

describe("routes", () => {
  it("round-trips a route and its parameters through the location hash", () => {
    const hash = formatRoute("trace", { group: "edge", name: "3<-1" });
    expect(hash).toBe("#/trace?group=edge&name=3%3C-1");
    const { route, params } = parseRoute(hash);
    expect(route).toBe("trace");
    expect(params.get("name")).toBe("3<-1");
  });

  it("reads an empty hash as the default route", () => {
    expect(parseRoute("").route).toBe("");
    expect(parseRoute("#/").route).toBe("");
    expect(formatRoute("convergence")).toBe("#/convergence");
  });
});
