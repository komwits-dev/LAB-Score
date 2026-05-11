#!/usr/bin/env bash
set -euo pipefail

# Replace personal absolute dbCAN paths in README with generic/project-local examples.
perl -0pi -e 's#/media/mecob/komwit/db/dbcan#/path/to/dbcan_database#g' README.md

echo "README.md path placeholders updated."
echo "Now review, then run:"
echo "  git add README.md"
echo "  git commit -m 'Use generic dbCAN database path in README'"
echo "  git push origin main"
