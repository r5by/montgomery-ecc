from mont.typing import GFElementType, GFType


class GFmodp(GFType):
    '''Simple Z/pZ field wrapper for quick testing purpose'''
    def __init__(self, p):
        # if not self._is_prime(p):
        #     raise ValueError("p must be a prime number.")
        self.modulus = p
        super().__init__(p)  # setup the group order

    def __call__(self, value):
        return GFmodpElement(value % self.modulus, self)

    def _is_prime(self, n):
        if n <= 1:
            return False
        if n <= 3:
            return True
        if n % 2 == 0 or n % 3 == 0:
            return False
        i = 5
        while i * i <= n:
            if n % i == 0 or n % (i + 2) == 0:
                return False
            i += 6
        return True


class GFmodpElement(GFElementType):
    def __eq__(self, other: 'GFElementType') -> bool:
        if isinstance(other, int):
            return self.value == other  # a convenient verification of its value only

        if isinstance(other, GFmodpElement):
            return self.value == other.value and self.field.modulus == other.field.modulus

        return NotImplemented

    def __init__(self, value: int, field: GFmodp):
        self.value = value
        self.field = field

    def __add__(self, other):
        if isinstance(other, GFmodpElement):
            # !! Here + must be included in (), b.c. the priority of + vs. % is not defined!!
            return GFmodpElement((self.value + other.value) % self.field.modulus, self.field)
        elif isinstance(other, int):
            return GFmodpElement((self.value + other) % self.field.modulus, self.field)
        return NotImplemented

    def __radd__(self, other: int):
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, GFmodpElement):
            return GFmodpElement((self.value - other.value) % self.field.modulus, self.field)
        return NotImplemented

    def __mul__(self, other):
        if isinstance(other, GFmodpElement):
            return GFmodpElement((self.value * other.value) % self.field.modulus, self.field)
        return NotImplemented

    def __truediv__(self, other):
        if isinstance(other, GFmodpElement) and other.value != 0:
            # Find multiplicative inverse using extended Euclidean algorithm
            inverse = pow(other.value, self.field.modulus - 2, self.field.modulus)
            return GFmodpElement((self.value * inverse) % self.field.modulus, self.field)
        elif isinstance(other, int):
            # Design choice: if we are dividing a '1', we should treat it as the identity in the Field;
            # o.w. treat it an 'as-is' integer already in Field repr.
            other = self.field(other) if other == 1 else GFmodpElement(other, self.field)
            return self.__truediv__(other)
        return NotImplemented

    def __rtruediv__(self, other: int):
        t = self.field(other) if other == 1 else GFmodpElement(other, self.field)
        return t.__truediv__(self)

    def __neg__(self):
        return GFmodpElement(-self.value % self.field.modulus, self.field)

    def __pow__(self, exponent):
        return GFmodpElement(pow(self.value, exponent, self.field.modulus), self.field)

    def __int__(self):
        return self.value

    def __repr__(self):
        return f"{self.value} mod {self.field.modulus}"


