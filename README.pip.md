# Install on LAM's cluster using pip

To install via pip, you must activate a virtual environment
```
python3.12 -m virtualenv ~/cesam/venv/amazed
source  ~/cesam/venv/amazed/bin/activate
```

## Debug mode

```
source start_amazed.sh
pip install -v . -C build-dir=build-pip -C cmake.build-type=Debug  -C cmake.define.CMAKE_PREFIX_PATH=$SHARED
```

## Debug mode, no TEST

```
source start_amazed.sh
pip install -v . -C build-dir=build-pip -C cmake.build-type=Debug -C cmake.define.BUILD_TESTING=OFF -C cmake.define.CMAKE_PREFIX_PATH=$SHARED
```

## Debug mode, editable ( -e ) and no TEST

```
source start_amazed.sh
pip install -v -e . -C build-dir=build-pip -C cmake.build-type=Debug -C cmake.define.BUILD_TESTING=OFF -C cmake.define.CMAKE_PREFIX_PATH=$SHARED
```

## Relase mode
```
source start_amazed.sh
pip install -v . -C build-dir=build-pip  -C cmake.define.CMAKE_PREFIX_PATH=$SHARED
```


# Test import module
You must be outside of your git directory. You can run :
```
cd
python -c "import pylibamazed.redshift as pr; print(pr.get_version()); import pylibamazed; print(pylibamazed.version)"
```

# Create a wheel
```
pip install build
python -m build -C build-dir=build-pip -C cmake.build-type=Debug -C cmake.define.BUILD_TESTING=OFF -C cmake.define.CMAKE_PREFIX_PATH=$SHARED
```
