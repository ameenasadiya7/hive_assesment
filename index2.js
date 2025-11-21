// Fully safe rational solver using BigInt fractions

const fs = require("fs");
const path = require("path");

const DIGITS = "0123456789abcdefghijklmnopqrstuvwxyz";

// Convert base-X string to BigInt
function parseInBase(s, base) {
    s = s.toLowerCase();
    let result = 0n;
    const B = BigInt(base);
    for (const ch of s) {
        const v = DIGITS.indexOf(ch);
        if (v < 0) throw new Error("Invalid digit: " + ch);
        result = result * B + BigInt(v);
    }
    return result;
}

// gcd for BigInt
function gcd(a, b) {
    while (b !== 0n) {
        let t = a % b;
        a = b;
        b = t;
    }
    return a < 0n ? -a : a;
}

// Rational number {num, den}
function makeFrac(num, den = 1n) {
    if (den === 0n) throw new Error("Zero denominator");
    let g = gcd(num, den);
    num /= g;
    den /= g;
    if (den < 0n) {
        num = -num;
        den = -den;
    }
    return { num, den };
}

function add(a, b) {
    return makeFrac(a.num * b.den + b.num * a.den, a.den * b.den);
}
function sub(a, b) {
    return makeFrac(a.num * b.den - b.num * a.den, a.den * b.den);
}
function mul(a, b) {
    return makeFrac(a.num * b.num, a.den * b.den);
}
function div(a, b) {
    return makeFrac(a.num * b.den, a.den * b.num);
}

// Gaussian elimination using Rational Arithmetic
function solve(matrix, rhs) {
    const n = matrix.length;

    // Forward elimination
    for (let col = 0; col < n; col++) {
        // Find pivot
        let pivot = col;
        for (let r = col; r < n; r++) {
            if (matrix[r][col].num !== 0n) {
                pivot = r; break;
            }
        }

        // Swap rows if needed
        if (pivot !== col) {
            [matrix[col], matrix[pivot]] = [matrix[pivot], matrix[col]];
            [rhs[col], rhs[pivot]] = [rhs[pivot], rhs[col]];
        }

        // Eliminate below
        for (let r = col + 1; r < n; r++) {
            if (matrix[r][col].num === 0n) continue;

            const factor = div(matrix[r][col], matrix[col][col]);

            for (let c = col; c < n; c++) {
                matrix[r][c] = sub(matrix[r][c], mul(factor, matrix[col][c]));
            }
            rhs[r] = sub(rhs[r], mul(factor, rhs[col]));
        }
    }

    // Back substitution
    const x = Array(n).fill(makeFrac(0n, 1n));
    for (let i = n - 1; i >= 0; i--) {
        let sum = makeFrac(0n);
        for (let j = i + 1; j < n; j++) {
            sum = add(sum, mul(matrix[i][j], x[j]));
        }
        x[i] = div(sub(rhs[i], sum), matrix[i][i]);
    }

    return x;
}

function main() {
    const data = JSON.parse(fs.readFileSync("./input.json", "utf8"));

    const k = data.keys.k;

    // Decode first k points
    const points = [];
    for (let i = 1; i <= k; i++) {
        const base = Number(data[i].base);
        const val = data[i].value;
        const y = parseInBase(val, base);
        points.push({ x: BigInt(i), y });
    }

    // Build Vandermonde matrix
    const A = Array.from({ length: k }, () =>
        Array(k).fill(makeFrac(0n))
    );
    const b = Array(k).fill(makeFrac(0n));

    for (let i = 0; i < k; i++) {
        let xpow = 1n;
        for (let j = 0; j < k; j++) {
            A[i][j] = makeFrac(xpow);
            xpow *= points[i].x;
        }
        b[i] = makeFrac(points[i].y);
    }

    // Solve A·coeff = b
    const coeff = solve(A, b);

    // constant term = coeff[0]
    const C = coeff[0];

    console.log(C.num / C.den); // should be an integer
}

main();
