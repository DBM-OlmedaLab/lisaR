# Risk-Based Test Suites

Suite layout establishes stable locations for focused tests. Production behavior is
not changed by this scaffold.

- `unit/`: isolated contract and validation tests.
- `integration/`: small synthetic pipeline-flow fixtures.
- `security/`: adversarial path and subprocess cases.
- `report-contract/`: semantic report and manifest assertions.
- `network-cache/`: offline cache fixtures; no live network in normal tests.
