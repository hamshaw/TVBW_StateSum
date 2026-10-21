jq -r '.order[]' Rank1-7.json |
    while IFS= read -r fruit; do
        python3 TVBW.py "$fruit" "Rank1-7.json" >> output_51.txt
    done
