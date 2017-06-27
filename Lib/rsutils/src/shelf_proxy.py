import shelve
from VVutil.Proxy import Proxy


class InvalidationError(Exception):
    pass

class ShelfProxy(Proxy):
    __slots__ = ["_key", "_shelf", "_invalidated"]
    
    def __init__(self, obj, shelf, key):
        Proxy.__init__(self, obj)
        object.__setattr__(self, "_shelf", shelf)
        object.__setattr__(self, "_key", key)
        object.__setattr__(self, "_invalidated", False)
    
    def __del__(self):
        try:
            sync_proxy(self)
        except InvalidationError:
            pass

class ShelfWrapper(object):
    def __init__(self, shelf):
        self.__shelf = shelf
        self.__cache = {}
    
    def __del__(self):
        self.close()
        
    def __getattr__(self, name):
        return getattr(self.__shelf, name)
    
    def __contains__(self, key):
        return key in self.__shelf
    
    def __len__(self, key):
        return len(self.__shelf)
    
    def __delitem__(self, key):
        if key in self.__cache:
            object.__setattr__(self.__cache[key], "_invalidated", True)
            del self.__cache[key]
        del self.__shelf[key]
    
    def __getitem__(self, key):
        try:
            obj = self.__cache[key]
        except KeyError:
            self.__cache[key] = obj = ShelfProxy(self.__shelf[key], self.__shelf, key)
        return obj
    
    def __setitem__(self, key, value):
        if key in self.__cache:
            object.__setattr__(self.__cache[key], "_invalidated", True)
            self.__cache[key] = ShelfProxy(value, self.__shelf, key)
        self.__shelf[key] = value
    
    def sync(self):
        for obj in self.__cache.itervalues():
            try:
                sync_proxy(obj)
            except InvalidationError:
                pass
    
    def close(self):
        self.sync()
        self.__cache.clear()
        self.__shelf.close()


def sync_proxy(proxy):
    if object.__getattribute__(proxy, "_invalidated"):
        raise InvalidationError("the proxy has been invalidated (the key was reassigned)")
    shelf = object.__getattribute__(proxy, "_shelf")
    key = object.__getattribute__(proxy, "_key")
    obj = object.__getattribute__(proxy, "_obj")
    shelf[key] = obj
    shelf.sync()

def open(*args):
    return ShelfWrapper( shelve.open(*args) )
