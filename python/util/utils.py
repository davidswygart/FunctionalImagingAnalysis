#TODO: should these always take ownership of the backing data?
# should we just use dataframes?

# # Provides the following notation:

# SoA[key: str] -> returns a list
# SoA[ind: int] -> returns a dict
# SoA[(key,ind) -> Tuple(str,int)] -> returns a scalar

class SoA:
    def __init__(self, dict_of_arrays):
        self._dict = dict_of_arrays
        self._keys = tuple(dict_of_arrays.keys())

    def keys(self):
        return self._keys

    def __len__(self):
        return len(self._dict[self._keys[0]])

    def __iter__(self):
        i = 0
        while i < len(self):
            yield self[i]
            i += 1

    def __getitem__(self, key):
        if isinstance(key, str):
            return self._dict[key]
        elif isinstance(key, int):
            return {k:self._dict[k][key] for k in self._keys}            
        elif isinstance(key, tuple):
            if isinstance(key[0], str) and isinstance(key[1], int):
                return self._dict[key[0]][key[1]]
            elif isinstance(key[1], str) and isinstance(key[0], int):
                return self._dict[key[1]][key[0]]
            else:
                raise Exception(f'tuple index must be of type (str,int) or (int,str), but was {tuple(map(type, key))}')
        else:
            raise Exception('index must be of type str, int, or tuple')

    
    def __setitem__(self, key, value):
        print(f'Key is type {type(key)}, value is type {type(value)}')
        raise NotImplementedError()

    def toAoS(self):
        return AoS(list(self.__iter__()))

    def copy(self):
        return SoA(self._dict.copy())

# Same as above, but backed by a different storage mechanism
class AoS:
    def __init__(self, array_of_dicts):
        self._array = array_of_dicts
        self._keys = tuple(self._array[0].keys())
        
    def keys(self):
        return self._keys
    
    def __len__(self):
        return len(self._array)
    
    def __iter__(self):
        for d in self._array:
            yield d
    
    def __getitem__(self, key):
        if isinstance(key, str):
            return [d[key] for d in self]
        elif isinstance(key, int):
            return self._array[key]        
        elif isinstance(key, tuple):
            if isinstance(key[0], str) and isinstance(key[1], int):
                return self._array[key[1]][key[0]]
            elif isinstance(key[1], str) and isinstance(key[0], int):
                return self._array[key[0]][key[1]]
            else:
                raise Exception(f'tuple index must be of type (str,int) or (int,str), but was {tuple(map(type, key))}')
        else:
            raise Exception('index must be of type str, int, or tuple')

    def __setitem__(self, key, value):
        raise NotImplementedError()

    def toSoA(self):
        return SoA({k:self[k] for k in self._keys})

    def copy(self):
        return SoA(self._array.copy())