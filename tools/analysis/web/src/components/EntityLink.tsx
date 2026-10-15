// An Infection, edge, pair, quantity or stratum as a link that selects it into the Inspector without leaving the
// workspace. It has a real href (the current URL with `sel`), so it can be opened in a new tab or copied.

import type { MouseEvent, ReactNode } from "react";
import { useLocation } from "@tanstack/react-router";
import { useSelection } from "../app/selection";
import { entityLabel, formatEntity, type Entity } from "../lib/entity";

export interface EntityLinkProps {
  entity: Entity;
  /** What the link shows; the entity's label by default. */
  children?: ReactNode;
  tab?: string;
  className?: string;
  title?: string;
}

export function EntityLink({ entity, children, tab, className = "", title }: EntityLinkProps) {
  const { entity: selected, select, searchFor } = useSelection();
  const { pathname } = useLocation();
  const active = selected !== null && formatEntity(selected) === formatEntity(entity);
  const mono = entity.kind === "infection" || entity.kind === "edge" || entity.kind === "pair";
  const onClick = (e: MouseEvent) => {
    if (e.metaKey || e.ctrlKey || e.shiftKey || e.button !== 0) return;
    e.preventDefault();
    e.stopPropagation();
    select(entity, tab);
  };
  return (
    <a
      href={pathname + searchFor(entity)}
      onClick={onClick}
      title={title ?? `Inspect ${entityLabel(entity)}`}
      data-entity={formatEntity(entity)}
      aria-current={active || undefined}
      className={`rounded-[2px] text-ink underline decoration-rule decoration-1 underline-offset-2 hover:text-accent hover:decoration-accent ${mono ? "id" : ""} ${active ? "bg-selected" : ""} ${className}`}
    >
      {children ?? entityLabel(entity)}
    </a>
  );
}
