if [[ `uname` == "Linux" ]]; then
    export LDSHARED="-shared -pthread"
else
    export LDSHARED="-bundle -undefined dynamic_lookup"
fi
$PYTHON -m pip install . --no-deps -vv
