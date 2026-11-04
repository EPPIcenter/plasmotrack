import { defineConfig } from "vitest/config";

// `pnpm dev` serves the frontend with hot reload and forwards /api to a running `pt-analysis serve` (port 8765).
export default defineConfig({
  server: { proxy: { "/api": "http://127.0.0.1:8765" } },
  build: { outDir: "dist", emptyOutDir: true },
  test: { include: ["src/**/*.test.ts"], environment: "node" },
});
