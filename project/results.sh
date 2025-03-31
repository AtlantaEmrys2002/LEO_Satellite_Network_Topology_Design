# USE TO GENERATE OVERVIEW OF RESULTS
echo "RESULTS OVERVIEW"
echo "Building table..."

# Format results
python ./analysis/build_results_table.py

echo "Table completed and stored in results_overview.csv"
