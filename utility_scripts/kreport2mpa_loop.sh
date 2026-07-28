for file in *report*; do
    # Create the output filename by replacing '.report' with '.mpa'
    outfile="${file%.report}.mpa"
    
    # Execute the python script
    python ../../utility_scripts/kreport2mpa.py -r "$file" -o "$outfile" --display-header
    
    echo "Converted $file -> $outfile"
done