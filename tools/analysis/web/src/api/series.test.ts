import { afterEach, describe, expect, it, vi } from "vitest";
import { failure, mockApi } from "../test/render";
import { loadSeries } from "./series";

const qs = (url: string) => new URL(url, "http://test").searchParams;

afterEach(() => vi.unstubAllGlobals());

describe("quantities' compact series", () => {
  it("go as one request when asked for within a few milliseconds, each answered with its own Chains' series", async () => {
    const fetch = mockApi();
    const [logPosterior, meanQ] = await Promise.all([
      loadSeries({ group: "log posterior", name: "log posterior" }),
      loadSeries({ group: "functionals", name: "mean_q" }),
    ]);
    expect(fetch).toHaveBeenCalledTimes(1);
    const params = qs(String(fetch.mock.calls[0]![0]));
    expect(params.getAll("q")).toEqual(["log posterior:log posterior", "functionals:mean_q"]);
    expect(params.get("points")).toBe("40");
    // fixtures/series.json: two Chains of 10 Draws each.
    expect(logPosterior?.map((c) => c.length)).toEqual([10, 10]);
    expect(meanQ).toHaveLength(2);
  });

  it("is null for a quantity no Chain has", async () => {
    mockApi();
    expect(await loadSeries({ group: "edge", name: "9<-1" })).toBeNull();
  });

  it("asks at most 200 quantities per request", async () => {
    const fetch = mockApi({ "/api/series": (url: string) => ({ chains: ["chain1"], points: 40, series: qs(url).getAll("q").map((q) => ({ group: "g", name: q.slice(2), chains: [[1]] })) }) });
    const all = await Promise.all(Array.from({ length: 250 }, (_, i) => loadSeries({ group: "g", name: String(i) })));
    expect(fetch.mock.calls.map(([u]) => qs(String(u)).getAll("q").length)).toEqual([200, 50]);
    expect(all.every((s) => s?.[0]?.[0] === 1)).toBe(true);
  });

  it("fails every quantity of a request that failed", async () => {
    mockApi({ "/api/series": failure(400, "a quantity is GROUP:NAME") });
    await expect(loadSeries({ group: "x", name: "y" })).rejects.toThrow("a quantity is GROUP:NAME");
  });
});
