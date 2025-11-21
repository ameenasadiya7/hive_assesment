
const fs = require('fs');
const path = require('path');

function parseInput(filePath) {
  const raw = fs.readFileSync(filePath, 'utf8');
  const obj = JSON.parse(raw);
  const n = obj.keys && obj.keys.n ? Number(obj.keys.n) : null;
  const k = obj.keys && obj.keys.k ? Number(obj.keys.k) : null;
  if (n === null || k === null) {
    throw new Error("JSON must contain keys.n and keys.k");
  }

  
  const points = [];
  Object.keys(obj).forEach(kname => {
    if (kname === 'keys') return;
    
    if (!/^\d+$/.test(kname)) return;
    const x = Number(kname);
    const entry = obj[kname];
    if (!entry || !entry.base || entry.value === undefined) {
      throw new Error(`Invalid entry for key ${kname}`);
    }
    
    const baseNum = Number(entry.base);
    if (!Number.isInteger(baseNum) || baseNum < 2 || baseNum > 36) {
      throw new Error(`Unsupported base ${entry.base} for key ${kname}`);
    }
    const y = parseInt(entry.value, baseNum);
    if (Number.isNaN(y)) {
      throw new Error(`Could not parse value '${entry.value}' in base ${entry.base} for key ${kname}`);
    }
    points.push({ x, y });
  });

  if (points.length < k) {
    throw new Error(`Not enough points provided (${points.length}) for required k = ${k}`);
  }

  
  points.sort((a, b) => a.x - b.x);

  return { n, k, points };
}


function buildSystem(points, degree) {
  const k = degree + 1;
  const A = new Array(k).fill(0).map(() => new Array(k).fill(0));
  const b = new Array(k).fill(0);
  for (let i = 0; i < k; i++) {
    const xi = points[i].x;
    let pow = 1;
    for (let j = 0; j < k; j++) {
      A[i][j] = pow;
      pow *= xi;
    }
    b[i] = points[i].y;
  }
  return { A, b };
}


function solveLinearSystem(A, b) {
  const n = A.length;
  
  const M = A.map(row => row.slice());
  const y = b.slice();

  for (let k = 0; k < n; k++) {
    
    let maxRow = k;
    let maxVal = Math.abs(M[k][k]);
    for (let r = k + 1; r < n; r++) {
      const val = Math.abs(M[r][k]);
      if (val > maxVal) {
        maxVal = val;
        maxRow = r;
      }
    }
    if (maxVal === 0) {
      throw new Error("Singular matrix encountered while solving linear system");
    }
    
    if (maxRow !== k) {
      [M[k], M[maxRow]] = [M[maxRow], M[k]];
      [y[k], y[maxRow]] = [y[maxRow], y[k]];
    }

    
    for (let i = k + 1; i < n; i++) {
      const factor = M[i][k] / M[k][k];
      for (let j = k; j < n; j++) {
        M[i][j] -= factor * M[k][j];
      }
      y[i] -= factor * y[k];
    }
  }

  
  const x = new Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    let sum = 0;
    for (let j = i + 1; j < n; j++) sum += M[i][j] * x[j];
    x[i] = (y[i] - sum) / M[i][i];
  }
  return x;
}

function main() {
  const inputFile = process.argv[2] || path.join(__dirname, 'input.json');
  try {
    const { n, k, points } = parseInput(inputFile);
    const degree = k - 1;

    
    const samplePoints = points.slice(0, k);
    const { A, b } = buildSystem(samplePoints, degree);

    const coeffs = solveLinearSystem(A, b); // coeffs[0] is constant term c0

    const C = coeffs[0];

    
    let ok = true;
    for (let i = k; i < points.length; i++) {
      const p = points[i];
      let val = 0;
      for (let j = 0; j <= degree; j++) {
        val += coeffs[j] * Math.pow(p.x, j);
      }
      if (Math.abs(val - p.y) > 1e-6) {
        ok = false;
        
      }
    }

    
    console.log(C);

    

  } catch (err) {
    console.error("Error:", err.message);
    process.exit(1);
  }
}

if (require.main === module) main();

