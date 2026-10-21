// A tooltip in the document body at a viewport point, so no panel, drawer or chart frame clips it: it measures itself
// and sits beside the point, flipped where it would cross the viewport's edge (lib/tip). Every chart's tip and the
// network's go through it.

import { useLayoutEffect, useRef, useState, type ReactNode } from "react";
import { createPortal } from "react-dom";
import { placeTip, type TipLine } from "../lib/tip";

export function FloatingTip({ point, children }: { point: { x: number; y: number }; children: ReactNode }) {
  const ref = useRef<HTMLDivElement>(null);
  const [at, setAt] = useState<{ left: number; top: number } | null>(null);

  useLayoutEffect(() => {
    const box = ref.current!.getBoundingClientRect();
    setAt(placeTip(point, { width: box.width, height: box.height }, { width: window.innerWidth, height: window.innerHeight }));
  }, [point.x, point.y, children]);

  return createPortal(
    <div
      ref={ref}
      role="tooltip"
      className="pointer-events-none z-50 min-w-36 max-w-96 rounded-[3px] border border-line bg-raised px-2 py-1 text-axis text-ink shadow-sm"
      style={{ position: "fixed", left: at?.left ?? point.x, top: at?.top ?? point.y, visibility: at ? "visible" : "hidden" }}
    >
      {children}
    </div>,
    document.body,
  );
}

/** A Plot tip's lines: a bold label and its value per line, the colour swatch after the value. */
export function TipLines({ lines }: { lines: readonly TipLine[] }) {
  return (
    <>
      {lines.map((line, i) => (
        <div key={i} className="whitespace-pre tabular-nums">
          {line.label && <span className="font-semibold">{line.label}</span>}
          {line.label && " "}
          {line.value}
          {line.swatch && (
            <span aria-hidden style={{ color: line.swatch }}>
              {" ■"}
            </span>
          )}
        </div>
      ))}
    </>
  );
}
