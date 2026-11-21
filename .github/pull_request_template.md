## What and why

<!-- What this changes, and the problem it solves. Link any issue. -->

## Output

- [ ] Output is unchanged (byte for byte against the previous build), or
- [ ] Output changes on purpose, as a bug fix listed in `CHANGELOG.md`
      under *Changes to output*

## Checks

- [ ] `make test` passes (and `make SANITIZE=1 test` for memory-related changes)
- [ ] New or changed behaviour has a unit test or a `tests/cli` fixture
- [ ] For a speed-up: timings before and after, with the machine
- [ ] Documentation updated (README.md, docs/)
