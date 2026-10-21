import { defineConfig, devices } from "@playwright/test";

// `pnpm e2e`: build, then serve the golden Analysis (scripts/golden.py) with the build and visit it in Chromium.
// Not part of ../test.sh.
const PORT = Number(process.env.PT_E2E_PORT ?? 8899);

export default defineConfig({
  testDir: "e2e",
  timeout: 30_000,
  reporter: "list",
  use: { baseURL: `http://127.0.0.1:${PORT}`, trace: "retain-on-failure" },
  projects: [{ name: "chromium", use: { ...devices["Desktop Chrome"], viewport: { width: 1440, height: 900 } } }],
  webServer: {
    command: `cd .. && uv run python web/scripts/golden.py serve --port ${PORT}`,
    url: `http://127.0.0.1:${PORT}/api/analysis`,
    reuseExistingServer: false,
    timeout: 120_000,
  },
});
