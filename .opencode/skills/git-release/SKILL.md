---
name: git-release
description: Create consistent releases and changelogs
---

## What I do
- Draft release notes from merged PRs
- Propose a version bump
- Provide a copy-pasteable `gh release create` command

## When to use me
Use this when preparing a tagged release.
Ask clarifying questions if versioning scheme is unclear.

## Process
1. Run `git log --oneline $(git describe --tags --abbrev=0)..HEAD`
2. Categorise commits (features, fixes, chores)
3. Determine version bump (major/minor/patch)
4. Generate changelog entry
5. Output the release command