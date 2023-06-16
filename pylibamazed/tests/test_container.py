from pylibamazed.Container import Container
import numpy as np

# Tests in case T in a numpy array
def test_container_equality():
    # Two empty containers are equal
    assert(Container() == Container())

    array1 = np.array([1,2,3])
    array2 = np.array([0,2,3])
    # Container 1 included in container 2: not equal
    assert(Container(**{"key": array1}) != Container(**{"key": array1, "key2": array1}))

    # Containers with same data but different key not equal
    assert(Container(**{"key": array1}) != Container(**{"other key": array1}))

    # Container with same keys but different data not equal
    assert(Container(**{"key": array1}) != Container(**{"other key": array2}))

    # Containers with same data and keys but in different order are not equal
    assert(Container(**{"key": array1, "key2": array2}) == Container(**{"key2": array2, "key": array1}))

