// A link to a workspace (or the Inspector's page) with its own parameters, keeping the selection, multi-selection and
// active Analysis (`linkTo`). It has a real href, so it opens in a new tab with a modifier key; a plain click
// navigates on the client router.

import type { MouseEvent, ReactNode } from "react";
import { useLinkTo, type Destination, type SearchValues } from "../app/search";

export interface WorkspaceLinkProps {
  to: Destination;
  search?: SearchValues;
  children: ReactNode;
  className?: string;
  title?: string;
  "aria-label"?: string;
  "aria-current"?: "page";
}

export function WorkspaceLink({ to, search, children, className = "text-accent hover:underline", title, "aria-label": label, "aria-current": current }: WorkspaceLinkProps) {
  const { href, go } = useLinkTo();
  const onClick = (e: MouseEvent) => {
    if (e.metaKey || e.ctrlKey || e.shiftKey || e.button !== 0) return;
    e.preventDefault();
    go(to, search);
  };
  return (
    <a href={href(to, search)} onClick={onClick} className={className} title={title} aria-label={label} aria-current={current}>
      {children}
    </a>
  );
}
