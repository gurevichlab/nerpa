from pathlib import Path
import polars as pl

def rust_f64_literal(mass_str: str) -> str:
    s = str(mass_str).strip()
    if s == '':
        raise ValueError('Empty AtomicMass value encountered')

    # Ensure Rust treats it as a float literal.
    if all(ch.isdigit() for ch in s):
        s = f'{s}.0'
    return s

csv_path = Path("/home/ilianolhin/work/tools/nerpa/tmp/atomic_masses.csv")
df = pl.read_csv(csv_path).select(['Symbol', 'AtomicMass'])

pairs: list[tuple[str, str]] = [(s, m) for s, m in df.iter_rows()]  # type: ignore[misc]

print("pub static atomic_masses: LazyLock<HashMap<&'static str, f64>> = LazyLock::new(|| {{")
print('    [')
for symbol, mass_str in pairs:
    lit = rust_f64_literal(mass_str)
    print(f'        ("{symbol}", {lit}),')
print('    ]')
print('    .into_iter()')
print('    .collect()')
print('});')

    
