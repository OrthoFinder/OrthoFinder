import functools
import time
import inspect

def timeit(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start = time.perf_counter()
        result = func(*args, **kwargs)
        end = time.perf_counter()
        time_elapsed = end - start
        print()
        print('***Function running time***')
        print(f"{func.__name__}: {time_elapsed:5f}s", end='\n'*2)
        return result
    return wrapper

