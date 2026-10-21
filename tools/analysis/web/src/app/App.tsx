// The providers: TanStack Query (cache keyed by endpoint and parameters), the theme, and the router.

import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { RouterProvider, type RouterHistory } from "@tanstack/react-router";
import { useState } from "react";
import { ApiError } from "../api/client";
import type { ThemeMode } from "../lib/palette";
import { createAppRouter } from "./router";
import { ThemeProvider } from "./theme";

/** A request the backend refused (a 4xx: unknown quantity, missing input) fails at once; others are tried twice. */
const retry = (failures: number, error: Error) => !(error instanceof ApiError && error.status < 500) && failures < 1;

export function createQueryClient() {
  return new QueryClient({ defaultOptions: { queries: { retry, refetchOnWindowFocus: false } } });
}

export function App({ history, theme }: { history?: RouterHistory; theme?: ThemeMode }) {
  const [client] = useState(createQueryClient);
  const [router] = useState(() => createAppRouter(history));
  return (
    <QueryClientProvider client={client}>
      <ThemeProvider mode={theme}>
        <RouterProvider router={router} />
      </ThemeProvider>
    </QueryClientProvider>
  );
}
