import sys
import json
import argparse


def searchForContext(context):
    print('Searching for a status with context: %s' % context)
    statuses = json.load(sys.stdin)

    for stat in statuses:
        if 'context' in stat and stat['context'] == context:
            sys.exit(0)

    sys.exit(1)


# =============================================================================
# Main
# =============================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser("Parse github api status list")
    parser.add_argument("--context", default=None, help="context of interest")
    args = parser.parse_args()
    searchForContext(args.context)
