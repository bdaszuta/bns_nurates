include config.mk

BIND_BUILD_DIR = bind_py/build

.PHONY: bind
bind:
	@test -n "$(PYTHON)" || (echo "PYTHON not set; activate conda first" && false)
	@$(PYTHON) -c "import nanobind" 2>/dev/null || \
		(echo "nanobind not found; install with: pip install nanobind" && false)
	cmake -S bind_py -B $(BIND_BUILD_DIR) \
		-DCMAKE_BUILD_TYPE=Release \
		-DPython_EXECUTABLE=$(PYTHON) \
		-Dnanobind_DIR=$(NANOBIND_CMAKE_DIR) \
		$(if $(GSL_ROOT),-DGSL_ROOT=$(GSL_ROOT),)
	cmake --build $(BIND_BUILD_DIR)

.PHONY: clean-bind
clean-bind:
	rm -rf $(BIND_BUILD_DIR)

.PHONY: test-bindings
test-bindings:
	bash tests/bind_py/run_all.sh
