#!/bin/bash
set -e -u -x

yum --disablerepo=epel -y update  ca-certificates

PLAT=manylinux2014_x86_64

function repair_wheel {
    wheel="$1"
    if ! auditwheel show "$wheel"; then
        echo "Skipping non-platform wheel $wheel"
    else
        auditwheel repair "$wheel" --plat "$PLAT" -w dist
    fi
}

rm -rf dist

# Install a system package required by our library
yum install -y atlas-devel

# Compile wheels
for PYBIN in /opt/python/*/bin; do
    "${PYBIN}/pip" install -r dev-requirements.txt
    "${PYBIN}/pip" wheel . --no-deps -w dist
done

# Bundle external shared libraries into the wheels
for whl in dist/*.whl; do
    repair_wheel "$whl"
done

/opt/python/cp36-cp36m/bin/pip install twine
/opt/python/cp36-cp36m/bin/twine upload --skip-existing dist/*manylinux*

