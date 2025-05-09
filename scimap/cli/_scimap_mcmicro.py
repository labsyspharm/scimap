# import scimap as sm
from .. import tools as tl
from .. import preprocessing as pp
from .. import plotting as pl
from .. import helpers as hl

import sys
import argparse
import pathlib
import textwrap
from joblib import Parallel, delayed
import time


# Actual mcmicro code


def mcmicro_wrap(argv=sys.argv):
    parser = argparse.ArgumentParser()

    parser.add_argument(
        'input_csv',
        help='single-cell quantification table (CSV file) generated by mcmicro',
    )
    parser.add_argument('-o', '--output', default='.', help='output directory')
    parser.add_argument(
        '--method',
        choices=['all', 'spatial', 'kmeans', 'phenograph', 'leiden'],
        default=['all'],
        nargs='+',
        help='specify one or more clustering methods to run',
    )
    parser.add_argument(
        '--csv',
        action='store_true',
        help='output csv version of the h5ad file',
    )

    args = parser.parse_args(argv[1:])

    # Run either all methods or specified methods
    methods = (
        ['spatial', 'kmeans', 'phenograph', 'leiden']
        if 'all' in args.method
        else args.method
    )

    results = []
    for m in methods:
        # Pass arguments as a list to clustering function
        result = clustering(
            [
                'clustering',  # Program name as first arg
                args.input_csv,  # The input CSV path
                '-o',
                args.output,
                '--clustering-method',
                m,
            ]
        )
        results.append(result)

    print(results)

    # Run merge if we have more than one method
    if len(methods) > 1:
        merge_args = [None, args.output, '-o', args.output]
        print("Running merge operation...")
        merge_result = merge(merge_args)  # Run merge first
        if merge_result != 0:
            print("Merge operation failed")
            return merge_result

        # Give filesystem time to sync after merge
        time.sleep(2)

        # move pdf files to output / plots
        output_dir = pathlib.Path(args.output)
        pdfs = sorted(output_dir.rglob('*.pdf'))
        plots_dir = output_dir / 'plots'
        plots_dir.mkdir(exist_ok=True)

        if len(pdfs) > 0:
            print(textwrap.indent('Moving pdf plots:\n', '    '))
        for p in pdfs:
            if plots_dir not in p.parents:
                move_to = plots_dir / p.parent.name / p.name
                print(textwrap.indent(f'{str(p)}\n    -> {str(move_to)}', '    '))
                move_to.parent.mkdir(exist_ok=True, parents=True)
                p.replace(move_to)
        print()

    # Generate CSV if requested, regardless of number of methods
    if args.csv:
        output_dir = pathlib.Path(args.output)
        try:
            csv_generated = False

            # For default/all methods case, check merged file first
            if 'all' in args.method or len(methods) > 1:
                merged_file = output_dir / 'combined_adata.h5ad'  # Updated filename
                print(f"\nLooking for merged file:")
                print(f"  Expected merged file: {merged_file}")
                print(f"  Directory contents:")
                for f in output_dir.glob('*'):
                    print(f"    {f}")
                if merged_file.exists():
                    print(
                        f"Found merged file, size: {merged_file.stat().st_size} bytes"
                    )
                    print(f"Generating CSV for merged file: {merged_file}")
                    try:
                        hl.scimap_to_csv(adata=str(merged_file), output_dir=output_dir)
                        csv_file = output_dir / f'{merged_file.stem}.csv'
                        if csv_file.exists():
                            print(f"Successfully created merged CSV: {csv_file}")
                            csv_generated = True
                        else:
                            print(f"Failed to create CSV for merged file")
                    except Exception as e:
                        print(f"Error generating CSV for merged file: {str(e)}")
                else:
                    print(f"Merged file not found at: {merged_file}")

            # Then check individual method directories
            for method in methods:
                method_dir = output_dir / method
                h5ad_files = list(method_dir.glob('*.h5ad'))
                for h5ad_file in h5ad_files:
                    if h5ad_file.exists():
                        print(f"Generating CSV for {h5ad_file}")
                        hl.scimap_to_csv(adata=str(h5ad_file), output_dir=method_dir)
                        csv_generated = True

            if not csv_generated:
                print("Warning: No h5ad files found to convert to CSV")

        except Exception as e:
            print(f"Error generating CSV: {str(e)}")

    return 0


def clustering(argv=sys.argv):
    parser = argparse.ArgumentParser()

    parser.add_argument(
        'csv',
        help='single-cell quantification table (CSV file) generated by mcmicro',
    )
    parser.add_argument('-o', '--output', default='.', help='output directory')
    clustering_methods = ['all', 'spatial', 'kmeans', 'phenograph', 'leiden', 'pass']
    parser.add_argument(
        '--clustering-method',
        default=['all'],
        choices=clustering_methods,
        nargs='+',
        help='choice of clustering algorithms, "pass": do not run any clustering methods; "all": run all clustering methods; default: "all"',
    )

    args = parser.parse_args(argv[1:])

    _output_dir = pathlib.Path(args.output)
    mcmicro_csv_path = args.csv
    methods = set(args.clustering_method)

    assert '.csv' in mcmicro_csv_path, 'input file must be a csv file'

    if 'all' in methods:
        methods = clustering_methods[1:-1]

    if 'pass' in methods:
        pp.mcmicro_to_scimap(
            feature_table_path=mcmicro_csv_path, output_dir=str(_output_dir)
        )
        return

    for method in methods:
        output_dir = _output_dir / method
        output_dir.mkdir(parents=True, exist_ok=True)  # Ensure directory exists

        print(f"Processing {method} clustering...")

        pp.mcmicro_to_scimap(
            feature_table_path=mcmicro_csv_path, output_dir=str(output_dir)
        )

        adata_path = (
            pathlib.Path(output_dir) / f'{pathlib.Path(mcmicro_csv_path).stem}.h5ad'
        )

        print(f"Created {adata_path}")

        if method == 'spatial':
            # Spatial clustering
            tl.spatial_expression(adata=str(adata_path), output_dir=output_dir)
            tl.spatial_cluster(
                adata=str(adata_path),
                df_name="spatial_expression",
                output_dir=output_dir,
            )

        else:
            # Expression clustering
            tl.cluster(adata=str(adata_path), method=method, output_dir=output_dir)

            # Expression clustering plots
            pl.cluster_plots(
                adata=str(adata_path), group_by=method, output_dir=output_dir
            )
    return 0


def merge(argv=sys.argv):
    parser = argparse.ArgumentParser()

    parser.add_argument(
        'directory',
        help='recursively search for .h5ad files to merge',
    )
    parser.add_argument('-o', '--output', default='.', help='output directory')
    parser.add_argument(
        '-d',
        '--delete-merged',
        default=False,
        action='store_true',
        help='delete found input files after merging; default: False',
    )
    parser.add_argument(
        '--csv',
        default=False,
        action='store_true',
        help='output csv version of the merged h5ad file; default: False',
    )

    args = parser.parse_args(argv[1:])
    print(f"Merge received args: {args}")  # Debug print

    input_dir = pathlib.Path(args.directory)
    output_dir = pathlib.Path(args.output)
    delete_after = args.delete_merged
    output_csv = args.csv

    print(f"output_csv flag value: {output_csv}")  # Debug print

    input_files = sorted(input_dir.rglob('*.h5ad'))
    print(f"\nFound input files:")
    for f in input_files:
        print(f"  {f}")

    output_file = output_dir / 'combined_adata.h5ad'
    print(f"\nExpected output file path: {output_file}")
    print(f"Output directory contents:")
    for f in output_dir.glob('*'):
        print(f"  {f}")

    if len(input_files) == 0:
        print(f'No .h5ad files found in {str(input_dir)}')
        return 1

    # Verify all files exist
    missing_files = [f for f in input_files if not f.exists()]
    if missing_files:
        print("Error: The following files are missing:")
        for f in missing_files:
            print(f"  {f}")
        return 1

    try:
        # Merge data
        print(f"Starting merge with csv flag: {output_csv}")

        # Create output directory if it doesn't exist
        output_dir.mkdir(parents=True, exist_ok=True)

        # Ensure all input files exist and are readable
        for f in input_files:
            if not f.exists():
                print(f"Error: Input file not found: {f}")
                return 1

        # Perform merge operation
        print("Merging h5ad files...")
        try:
            # Store the result from merge_adata_obs
            merged_adata = hl.merge_adata_obs(
                adata=[str(p) for p in input_files], output_dir=output_dir
            )
            print(f"Merge operation completed. Result type: {type(merged_adata)}")

            # Update output file path to use correct filename
            output_file = output_dir / 'combined_adata.h5ad'

            print("\nChecking for merged file:")
            print(f"  Looking for: {output_file}")
            print(f"  Exists: {output_file.exists()}")

            print("\nFull directory contents after merge:")
            for f in output_dir.rglob('*'):
                print(f"  {f}")

        except Exception as merge_error:
            print(f"Error during merge_adata_obs: {str(merge_error)}")
            import traceback

            traceback.print_exc()

        # Wait and check for file with retries
        max_retries = 5
        retry_delay = 5  # seconds
        for attempt in range(max_retries):
            if output_file.exists():
                print(f"Merged file created successfully: {output_file}")
                break
            else:
                if attempt < max_retries - 1:
                    print(
                        f"Waiting for file to appear (attempt {attempt + 1}/{max_retries})..."
                    )
                    time.sleep(retry_delay)
                else:
                    print(
                        f"Error: Failed to create merged file after {max_retries} attempts"
                    )
                    return 1

        # Generate CSV if requested
        if output_csv:
            print(f"Generating CSV file for: {output_file}")
            hl.scimap_to_csv(adata=str(output_file), output_dir=output_dir)
            csv_file = output_file.parent / f'{output_file.stem}.csv'
            if csv_file.exists():
                print(f"Successfully created CSV file: {csv_file}")
            else:
                print(f"Warning: CSV file was not created: {csv_file}")

        if delete_after:
            if output_file in input_files:
                idx = input_files.index(output_file)
                input_files.pop(idx)
            print(
                textwrap.indent(
                    'Deleting:\n\n' + "\n".join([str(f) for f in input_files]), '    '
                )
            )
            for f in input_files:
                f.unlink()
            print()

        return 0  # Success

    except Exception as e:
        print(f"Error during merge: {str(e)}")
        import traceback

        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(mcmicro_wrap())
