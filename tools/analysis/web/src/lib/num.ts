// Numbers from the backend's JSON.

/** The backend writes non-finite numbers as the strings "Infinity", "-Infinity" and "NaN"; this reads them back. */
export function num(x: number | string | null): number | null {
  return typeof x === "string" ? Number(x) : x;
}
