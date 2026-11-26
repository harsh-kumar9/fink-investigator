PYTHON ?= python3
.PHONY: test replay verify verify-compass

test:
	$(PYTHON) -m unittest discover -s tests -v

replay:
	$(PYTHON) pilot.py replay --out results

verify-compass:
	$(PYTHON) experiments/compass-2025/verify.py
	$(PYTHON) experiments/compass-2025/batch-2/verify.py
	$(PYTHON) experiments/compass-2025/batch-3/verify.py
	$(PYTHON) experiments/compass-2025/batch-4/verify.py

verify: test
	@set -eu; tmp=$$(mktemp -d); trap 'rm -r "$$tmp"' EXIT; \
	$(PYTHON) pilot.py replay --out "$$tmp"; \
	for file in "$$tmp"/*; do diff -r "$$file" "results/$$(basename "$$file")"; done; \
	$(PYTHON) agents/review.py plan --packet "$$tmp/review_packet.json"
