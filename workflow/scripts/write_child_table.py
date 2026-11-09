"""Write the child-to-family lookup consumed by the STR summary."""
import csv
import json
import sys


def main():
    mapping = json.loads(sys.argv[1])
    with open(sys.argv[2], "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["child", "trio_id"])
        writer.writeheader()
        for family, child in mapping.items():
            writer.writerow({"trio_id": family, "child": child})


if __name__ == "__main__":
    main()
