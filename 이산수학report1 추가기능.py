def matrixout(m, size):
    print("+" + "-" * (size * 11) + "+")
    for i in range(size):
        print("|", end="")
        for j in range(size):
            print("%10.3f" % m[i][j], end=" ")
        print("|")
    print("+" + "-" * (size * 11) + "+")

def getMatrixTranspose(m):
    return [[m[j][i] for j in range(len(m))] for i in range(len(m))]

def getMatrixMinor(m, i, j):
    return [
        row[:j] + row[j+1:]
        for r, row in enumerate(m) if r != i
    ]

def getMatrixDeterminant(m):
    n = len(m)

    if n == 0:
        return 1

    if n == 1:
        return m[0][0]

    if n == 2:
        return m[0][0] * m[1][1] - m[0][1] * m[1][0]

    determinant = 0

    for c in range(n):
        minor = getMatrixMinor(m, 0, c)
        determinant += (
            (-1) ** c
            * m[0][c]
            * getMatrixDeterminant(minor)
        )

    return determinant

def has_inverse(m):
    determinant = getMatrixDeterminant(m)

    if abs(determinant) < 1e-10:
        return False

    return True

def getMatrixInverse(m):
    n = len(m)
    determinant = getMatrixDeterminant(m)

    if abs(determinant) < 1e-10:
        return None

    cofactors = []

    for r in range(n):
        cofactorRow = []

        for c in range(n):
            minor = getMatrixMinor(m, r, c)
            cofactor = (
                (-1) ** (r + c)
                * getMatrixDeterminant(minor)
            )
            cofactorRow.append(cofactor)

        cofactors.append(cofactorRow)

    adjugate = getMatrixTranspose(cofactors)

    inverse = [
        [adjugate[r][c] / determinant for c in range(n)]
        for r in range(n)
    ]

    return inverse

def getMatrixInverseGaussJordan(m):
    n = len(m)

    a = [
        [float(x) for x in m[i]]
        + [float(i == j) for j in range(n)]
        for i in range(n)
    ]

    for i in range(n):
        pivot = max(range(i, n), key=lambda r: abs(a[r][i]))

        if abs(a[pivot][i]) < 1e-10:
            return None

        if pivot != i:
            a[i], a[pivot] = a[pivot], a[i]

        pivot_value = a[i][i]
        a[i] = [x / pivot_value for x in a[i]]

        for r in range(n):
            if r != i:
                factor = a[r][i]
                a[r] = [
                    a[r][c] - factor * a[i][c]
                    for c in range(2 * n)
                ]

    return [row[n:] for row in a]

def compareMatrix(a, b, tolerance=1e-8):
    if a is None or b is None:
        return False

    n = len(a)

    for i in range(n):
        for j in range(n):
            if abs(a[i][j] - b[i][j]) > tolerance:
                return False

    return True

def solveLinearSystem(inverse, b):
    n = len(inverse)
    x = [0.0] * n

    for i in range(n):
        for j in range(n):
            x[i] += inverse[i][j] * b[j]

    return x


def main():
    try:
        n = int(input("정방행렬의 차수를 입력하세요: "))

        if n <= 0:
            print("차수는 양의 정수여야 합니다.")
            return

        matrix = []

        for i in range(n):
            while True:
                try:
                    row = input(
                        str(i + 1) + "행의 원소를 입력하세요: "
                    ).split()

                    if len(row) != n:
                        print(
                            "원소를 정확히 "
                            + str(n)
                            + "개 입력하세요."
                        )
                        continue

                    matrix.append([float(x) for x in row])
                    break

                except ValueError:
                    print("숫자만 입력하세요.")

        print("\n입력한 행렬:")
        matrixout(matrix, n)

        determinant = getMatrixDeterminant(matrix)
        print("\n행렬식:", determinant)

        if not has_inverse(matrix):
            print("역행렬이 존재하지 않습니다.")
            return

        inverse1 = getMatrixInverse(matrix)
        inverse2 = getMatrixInverseGaussJordan(matrix)

        print("\n방법1 - 행렬식을 이용한 역행렬:")
        matrixout(inverse1, n)

        print("\n방법2 - 가우스-조던 소거법을 이용한 역행렬:")
        matrixout(inverse2, n)

        print("\n[결과 비교]")
        if compareMatrix(inverse1, inverse2):
            print("두 방법의 결과가 동일합니다.")
        else:
            print("두 방법의 결과가 다릅니다.")

        print("\n[연립방정식 자동 풀이]")
        print("상수항 벡터 b의 원소를 입력하세요.")

        while True:
            try:
                b = list(map(float, input(
                    "b의 원소를 입력하세요: "
                ).split()))

                if len(b) != n:
                    print(
                        "원소를 정확히 "
                        + str(n)
                        + "개 입력하세요."
                    )
                    continue

                break

            except ValueError:
                print("숫자만 입력하세요.")

        x = solveLinearSystem(inverse1, b)

        print("\n연립방정식의 해:")
        for i in range(n):
            print("x" + str(i + 1) + " =", round(x[i], 6))

    except ValueError:
        print("차수는 정수로 입력해야 합니다.")


if __name__ == "__main__":
    main()