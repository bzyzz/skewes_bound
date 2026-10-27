PYTHON   ?= python
WORKERS  ?= 48
CERT     = certification
RACE     = race
RACE_TS  = $(shell seq -s, 5000 5000 100000)

.PHONY: all certify s_cert smith check check-only arb-recert arb-combine bounds test clean clean-cache \
        race-build race-check race-zeros race-certify team-check

all: certify

certify: check

s_cert:
	cd $(CERT) && $(PYTHON) -u certify_S.py $(WORKERS)

smith:
	cd $(CERT) && $(PYTHON) -u certify_smith_unconditional.py

check: s_cert smith
	cd $(CERT) && $(PYTHON) -u check_theorem.py

# Re-check every shipped certificate against the Arb error bounds (no zero data needed).
check-only:
	cd $(CERT) && $(PYTHON) -u check_theorem.py
	cd $(CERT) && for c in s_cert_T5e7_arb s_certificate_arb s_cert_T8e7_st_arb; do $(PYTHON) -u check_theorem.py $$c.json || exit 1; done
	cd $(CERT) && $(PYTHON) -u certify_smith_unconditional.py

# Arb isolation of all zeros up to 8e7 (~15 h on 64 cores; needs the Arb zero caches, see README) and combination.
arb-recert:
	cd $(CERT) && $(PYTHON) -u arb_recert.py $(WORKERS)

arb-combine:
	cd $(CERT) && $(PYTHON) -u arb_combine.py
	cd $(CERT) && for c in s_cert_T5e7_arb s_certificate_arb s_cert_T8e7_st_arb; do $(PYTHON) -u check_theorem.py $$c.json || exit 1; done

bounds:
	cd $(CERT) && $(PYTHON) -u error_bounds.py

test:
	cd $(CERT) && $(PYTHON) -u test_certio.py

# ---- prime number race pi(x;24,11) - pi(x;24,1) ----
race-build:
	cd $(RACE)/criterion && $(CC) -O3 -fopenmp -shared -fPIC scan.c -o libscan.so -lm
	cd $(RACE)/heuristics && $(CC) -O3 -fopenmp -shared -fPIC scan.c -o libscan.so -lm

# Turing completeness of the shipped certified zeros, then the Arb certificate.
race-check:
	cd $(RACE) && for c in 3_2 4_3 8_5 24_11; do $(PYTHON) -u zeros/turing.py data/arb_$$c/certified.npz 40 $(RACE_TS) || exit 1; done
	$(MAKE) race-certify

race-certify:
	cd $(RACE)/criterion && $(PYTHON) -u certify.py 100000 111.5458300 1.2e-4 7.62e-5 40 0,0,1 $(WORKERS)

# Recompute the race zeros from scratch (float location, Arb certification, Turing; ~80 core-hours).
race-zeros:
	cd $(RACE) && PYTHON=$(PYTHON) zeros/run_float.sh > data/run_float.log 2>&1
	cd $(RACE) && PYTHON=$(PYTHON) zeros/run_certify.sh > data/run_certify.log 2>&1

clean:
	find . -name __pycache__ -type d -prune -exec rm -rf {} +

clean-cache:
	rm -rf $(CERT)/s_chunks

# Team races mod 8 and 12: Turing for the extra zeros (8.3, 12.11), then the Arb certificates.
team-check:
	cd $(RACE) && $(PYTHON) -u zeros/turing.py data/arb_8_3/certified.npz 40 $(RACE_TS)
	cd $(RACE) && $(PYTHON) -u zeros/turing.py data/arb_12_11/certified.npz 40 $(shell seq -s, 5000 5000 50000)
	cd $(RACE)/team && $(PYTHON) -u team_certify.py 8 100000 64.764612 1.4e-4 2e-5 40 0,0,1 $(WORKERS)
	cd $(RACE)/team && $(PYTHON) -u team_certify.py 12 50000 339.632308 2.8e-4 1.1e-4 60 0,0,1 $(WORKERS)
