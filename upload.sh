#!/usr/bin/env bash
# Automatic upload: stage -> commit -> rebase onto remote -> push
#
#   ./upload.sh                  commit message defaults to "auto update"
#   ./upload.sh "your message"
#
# Deliberately NOT using `git push -f`: a force push silently discards
# commits made from another machine (the server, another laptop).
set -euo pipefail

MSG="${1:-auto update}"
BRANCH="$(git rev-parse --abbrev-ref HEAD)"

echo "======== branch: $BRANCH ========"
git add -A
git status --short
echo

if git diff --cached --quiet; then
  echo "Nothing staged - working tree is clean."
else
  git commit -m "$MSG"
fi

echo
echo "======== syncing with origin/$BRANCH ========"
if ! git pull --rebase origin "$BRANCH"; then
  echo
  echo "Rebase stopped - most likely a conflict."
  echo "  fix the files, then:  git add <file> && git rebase --continue"
  echo "  to undo everything:   git rebase --abort"
  exit 1
fi

echo
echo "======== pushing ========"
git push origin "$BRANCH"
echo
echo "======== done ========"
