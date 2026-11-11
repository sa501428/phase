// Share live cutoff behavior across viewers without requesting the BAM again.
export function watchVaf(
  input,
  { current, busy = () => false, update, error },
) {
  let timer;
  const apply = () => {
    clearTimeout(timer);
    if (input.value === "" || !input.checkValidity()) return;
    const previous = current();
    if (previous == null) return;
    if (busy()) {
      timer = setTimeout(apply, 100);
      return;
    }
    const value = Number(input.value) / 100;
    if (value === previous) return;
    try {
      update(value);
    } catch (e) {
      error?.(e);
    }
  };
  input.addEventListener("input", () => {
    clearTimeout(timer);
    timer = setTimeout(apply, 200);
  });
  input.addEventListener("change", apply);
}
