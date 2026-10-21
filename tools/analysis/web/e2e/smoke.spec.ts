// Smoke test: every workspace route of the golden Analysis renders with no console errors, in both themes, and a
// selection opens the Inspector.

import { expect, test, type Page } from "@playwright/test";
import { workspaces } from "../src/app/workspaces";

function watchConsole(page: Page): string[] {
  const errors: string[] = [];
  page.on("console", (m) => m.type() === "error" && errors.push(m.text()));
  page.on("pageerror", (e) => errors.push(e.message));
  return errors;
}

for (const theme of ["light", "dark"] as const) {
  test(`every workspace renders without console errors (${theme})`, async ({ page }) => {
    const errors = watchConsole(page);
    await page.addInitScript((t) => localStorage.setItem("pt-analysis:theme", t), theme);
    for (const w of workspaces) {
      await page.goto(`/${w.id}`);
      await expect(page.getByRole("heading", { name: "golden" })).toBeVisible();
      await expect(page.getByRole("main", { name: "Workspace" })).not.toHaveText(/^Loading…$/);
      await page.waitForLoadState("networkidle");
      expect(errors, `console errors on /${w.id}`).toEqual([]);
    }
  });
}

test("the Overview, Convergence and the Inspector work on reload", async ({ page }) => {
  const errors = watchConsole(page);
  await page.goto("/overview");
  await expect(page.getByText("Worst R-hat")).toBeVisible();
  await page.getByRole("region", { name: "Attention" }).getByRole("link", { name: "1 → 3" }).click();
  await expect(page).toHaveURL(/sel=edge:1%3E3/);
  await page.goto("/convergence?sel=quantity:log%20posterior/log%20posterior");
  const inspector = page.getByRole("complementary", { name: /Quantity log posterior/ });
  await expect(inspector.getByText("Pooled")).toBeVisible();
  await page.goto("/inspect?sel=infection:3");
  await expect(page.getByText("Expected observed parents")).toBeVisible();
  await page.waitForLoadState("networkidle");
  expect(errors).toEqual([]);
});
