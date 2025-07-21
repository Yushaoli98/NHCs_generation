import pandas as pd, os
from nhc_generator.config import settings
from nhc_generator.process import process_row
from concurrent.futures import ProcessPoolExecutor, as_completed

def main():
    os.makedirs(settings.OUTPUT_DIR, exist_ok=True)
    os.makedirs(settings.FAILED_DIR, exist_ok=True)
    df = pd.read_csv(settings.CSV_FILE).dropna(subset=['ID', 'Original_SMILES'])
    failed = []
    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(process_row, row._asdict()) for row in df.itertuples(index=False)]
        for i, f in enumerate(as_completed(futures), 1):
            try:
                f.result()
            except Exception as e:
                print(f"Failed at {i}: {e}")
                failed.append(row.ID)
    if failed:
        pd.DataFrame(failed).to_csv(settings.FAILED_CSV)

if __name__ == '__main__':
    main()
