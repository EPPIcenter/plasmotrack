import tailwindcss from "@tailwindcss/vite";
import react from "@vitejs/plugin-react";
import { defineConfig } from "vitest/config";

// `pnpm dev` serves the app with hot reload and forwards /api to a running `pt-analysis serve`: port 8765, or the
// URL in PT_ANALYSIS_API (e.g. `PT_ANALYSIS_API=http://127.0.0.1:8820 pnpm dev`).
export default defineConfig({
  plugins: [react(), tailwindcss()],
  server: { proxy: { "/api": process.env.PT_ANALYSIS_API ?? "http://127.0.0.1:8765" } },
  build: { outDir: "dist", emptyOutDir: true, chunkSizeWarningLimit: 1500 },
  test: {
    include: ["src/**/*.test.{ts,tsx}"],
    environment: "jsdom",
    setupFiles: ["src/test/setup.ts"],
  },
});
