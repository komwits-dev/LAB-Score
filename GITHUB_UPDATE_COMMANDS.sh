#!/usr/bin/env bash
set -euo pipefail

# Run this inside your local LAB-Score repository folder:
# cd /media/mecob/komwit/Lacto_project/LAB_SCORE_v1/Release/lab_score_v35

git add README.md LICENSE .gitignore install.sh docs/images/
git commit -m "Improve README, installation guide, and report screenshots"
git push origin main
