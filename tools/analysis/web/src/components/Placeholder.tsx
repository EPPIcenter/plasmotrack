// What a workspace or Inspector tab that is not built yet shows.

export function Placeholder({ title, description }: { title: string; description: string }) {
  return (
    <div className="flex h-full min-h-40 items-center justify-center p-6" data-placeholder={title}>
      <div className="max-w-sm text-center">
        <p className="font-semibold">{title}</p>
        <p className="mt-1 text-muted">{description}</p>
        <p className="mt-2 text-axis text-muted">Not built yet in the React app.</p>
      </div>
    </div>
  );
}
