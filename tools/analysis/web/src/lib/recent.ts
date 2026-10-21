// Recent selections, newest first, for the command palette.

export const RECENT_LIMIT = 8;

/** `key` first, then the rest of `list` without it, at most `limit` entries. */
export function pushRecent(list: readonly string[], key: string, limit = RECENT_LIMIT): string[] {
  return [key, ...list.filter((k) => k !== key)].slice(0, limit);
}
