// Light and dark themes, both first class. The choice is saved; without one the system's is used (index.html sets
// it before the first paint).

import { createContext, useCallback, useContext, useLayoutEffect, useMemo, useState, type ReactNode } from "react";
import { THEMES, type ChartTheme, type ThemeMode } from "../lib/palette";

const KEY = "pt-analysis:theme";

interface ThemeValue {
  mode: ThemeMode;
  chart: ChartTheme;
  toggle: () => void;
}

const ThemeContext = createContext<ThemeValue>({ mode: "light", chart: THEMES.light, toggle: () => {} });

function initial(): ThemeMode {
  const set = document.documentElement.dataset.theme;
  if (set === "light" || set === "dark") return set;
  return window.matchMedia?.("(prefers-color-scheme: dark)").matches ? "dark" : "light";
}

export function ThemeProvider({ children, mode: fixed }: { children: ReactNode; mode?: ThemeMode }) {
  const [mode, setMode] = useState<ThemeMode>(() => fixed ?? initial());
  useLayoutEffect(() => {
    document.documentElement.dataset.theme = mode;
  }, [mode]);
  const toggle = useCallback(() => {
    setMode((m) => {
      const next = m === "light" ? "dark" : "light";
      try {
        localStorage.setItem(KEY, next);
      } catch {
        // private browsing: the choice lasts the session
      }
      return next;
    });
  }, []);
  const value = useMemo(() => ({ mode, chart: THEMES[mode], toggle }), [mode, toggle]);
  return <ThemeContext.Provider value={value}>{children}</ThemeContext.Provider>;
}

export const useTheme = () => useContext(ThemeContext);
