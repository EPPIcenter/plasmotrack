// What the Inspector shows: an Infection, an edge, a pair, a quantity or a stratum, and its text form in the URL's
// `sel` parameter.
//
//   infection:<id>                      edge:<parent>><child>            pair:<a>,<b>
//   quantity:<group>/<name>[@<stratum>] stratum:<column>[=<value>]
//
// A quantity's stratum (a functional of one epi stratum) is `@column=value`, `@column` for the Infections with no
// value of the column, or `@` alone for every Infection. A stratum without `=value` is the Infections with no value.

import type { Params } from "./params";

/** A functional's stratum: `column` null for every Infection; `value` null for the Infections with no value. */
export interface QuantityStratum {
  column: string | null;
  value: string | null;
}

export type Entity =
  | { kind: "infection"; id: string }
  | { kind: "edge"; parent: string; child: string }
  | { kind: "pair"; a: string; b: string }
  | { kind: "quantity"; group: string; name: string; stratum?: QuantityStratum }
  | { kind: "stratum"; column: string; value: string | null };

export type EntityKind = Entity["kind"];

/** Split at the first `sep`: [before, after], or null when `sep` is absent. */
function cut(text: string, sep: string): [string, string] | null {
  const i = text.indexOf(sep);
  return i < 0 ? null : [text.slice(0, i), text.slice(i + sep.length)];
}

function columnValue(text: string): { column: string; value: string | null } | null {
  const parts = cut(text, "=");
  const column = parts ? parts[0] : text;
  return column ? { column, value: parts ? parts[1] : null } : null;
}

/** The entity a `sel` value names, or null when it names none. */
export function parseEntity(text: string): Entity | null {
  const head = cut(text, ":");
  if (!head) return null;
  const [kind, rest] = head;
  switch (kind) {
    case "infection":
      return rest ? { kind, id: rest } : null;
    case "edge": {
      const p = cut(rest, ">");
      return p && p[0] && p[1] ? { kind, parent: p[0], child: p[1] } : null;
    }
    case "pair": {
      const p = cut(rest, ",");
      return p && p[0] && p[1] && p[0] !== p[1] ? { kind, a: p[0], b: p[1] } : null;
    }
    case "quantity": {
      const at = rest.lastIndexOf("@");
      const body = at < 0 ? rest : rest.slice(0, at);
      const p = cut(body, "/");
      if (!p || !p[0] || !p[1]) return null;
      const q: Entity = { kind, group: p[0], name: p[1] };
      if (at < 0) return q;
      const s = rest.slice(at + 1);
      return { ...q, stratum: s ? columnValue(s)! : { column: null, value: null } };
    }
    case "stratum": {
      const s = columnValue(rest);
      return s ? { kind, ...s } : null;
    }
    default:
      return null;
  }
}

/** The `sel` value of an entity; `parseEntity` reads it back. Also its identity, for comparing selections. */
export function formatEntity(e: Entity): string {
  switch (e.kind) {
    case "infection":
      return `infection:${e.id}`;
    case "edge":
      return `edge:${e.parent}>${e.child}`;
    case "pair":
      return `pair:${e.a},${e.b}`;
    case "quantity": {
      const s = e.stratum;
      const tail = !s ? "" : s.column === null ? "@" : `@${s.column}${s.value === null ? "" : `=${s.value}`}`;
      return `quantity:${e.group}/${e.name}${tail}`;
    }
    case "stratum":
      return `stratum:${e.column}${e.value === null ? "" : `=${e.value}`}`;
  }
}

function stratumText(s: QuantityStratum): string {
  if (s.column === null) return "all Infections";
  return `${s.column}: ${s.value ?? "no value"}`;
}

/** How the Inspector's header, links and the command palette name an entity. */
export function entityLabel(e: Entity): string {
  switch (e.kind) {
    case "infection":
      return e.id;
    case "edge":
      return `${e.parent} → ${e.child}`;
    case "pair":
      return `${e.a}, ${e.b}`;
    case "quantity":
      return `${e.group}: ${e.name}${e.stratum ? ` (${stratumText(e.stratum)})` : ""}`;
    case "stratum":
      return `${e.column} = ${e.value ?? "no value"}`;
  }
}

/** A quantity's parameters for /api/trace and /api/range: group, name, and its stratum as `stratum_column` (empty for
 *  every Infection) and `stratum` (left out for the Infections with no value). */
export function quantityParams(e: Extract<Entity, { kind: "quantity" }>): Params {
  const params: Params = { group: e.group, name: e.name };
  if (e.stratum) {
    params.stratum_column = e.stratum.column ?? "";
    if (e.stratum.column !== null && e.stratum.value !== null) params.stratum = e.stratum.value;
  }
  return params;
}
