(() => {
  'use strict';

  const DB = window.NN_PARAMETER_DB;
  let lastRows = [];
  const $ = id => document.getElementById(id);

  const complement = { A: 'U', U: 'A', G: 'C', C: 'G' };
  const suffix = { OMe: 'o', F: 'f', MOE: 'm', LNA: 'l' };
  // Internal modified-base symbols in the fitted parameter tables:
  // H = A, J = C, K = G, L = U.
  const internal = { A: 'H', C: 'J', G: 'K', U: 'L' };

  function tokenize(sequence) {
    const s = String(sequence || '')
      .replace(/5[′']/gi, '')
      .replace(/3[′']/gi, '')
      .replace(/[\s\-–—>]/g, '');

    const output = [];
    for (let i = 0; i < s.length;) {
      if ((s[i] === 'd' || s[i] === 'D') && i + 1 < s.length && /[ACGT]/i.test(s[i + 1])) {
        output.push('d' + s[i + 1].toUpperCase());
        i += 2;
        continue;
      }

      const base = s[i].toUpperCase();
      if (!/[ACGU]/.test(base)) {
        throw new Error(`Unsupported sequence symbol near “${s.slice(i, i + 4)}”.`);
      }

      if (i + 1 < s.length && /[ofml]/i.test(s[i + 1])) {
        output.push(base + s[i + 1].toLowerCase());
        i += 2;
      } else {
        output.push(base);
        i += 1;
      }
    }
    return output;
  }

  function canonicalBase(token) {
    if (token.startsWith('d')) return token === 'dT' ? 'U' : token[1];
    return token[0];
  }

  function isExpectedStrand2Token(mode, token) {
    if (mode === 'unmodified') return /^[ACGU]$/.test(token);
    if (mode === 'RNA_DNA') return /^d[ACGT]$/.test(token);
    return token.endsWith(suffix[mode]);
  }

  function strand2Token(mode, rnaBase) {
    if (mode === 'unmodified') return rnaBase;
    if (mode === 'RNA_DNA') return 'd' + (rnaBase === 'U' ? 'T' : rnaBase);
    return rnaBase + suffix[mode];
  }

  function generateComplement(showMessage = true) {
    try {
      const a = tokenize($('a').value);
      if (!a.length) throw new Error('Enter Strand 1 first.');
      if (!a.every(t => /^[ACGU]$/.test(t))) {
        throw new Error('Strand 1 must contain unmodified RNA residues A, C, G and U only.');
      }
      const mode = $('mode').value;
      $('b').value = a.map(t => strand2Token(mode, complement[canonicalBase(t)])).join('');
      if (showMessage) setMessage('Complementary strand generated.', true);
      else clearMessage();
    } catch (error) {
      setMessage(error.message, false);
    }
  }

  function makeStepKey(mode, topPair, bottomPair) {
    const top = topPair.map(canonicalBase).join('');
    if (mode === 'unmodified') {
      return `${top}/${bottomPair.map(canonicalBase).join('')}`;
    }
    if (mode === 'RNA_DNA') {
      const bottom = bottomPair.slice().reverse().map(t => canonicalBase(t) === 'U' ? 'T' : canonicalBase(t)).join('');
      return `r${top}/d${bottom}`;
    }
    return `${top}/${bottomPair.map(t => internal[canonicalBase(t)]).join('')}`;
  }

  function reverseString(s) {
    return s.split('').reverse().join('');
  }

  function isSelfComplementary(tokens) {
    const sequence = tokens.map(canonicalBase);
    const reverseComplement = sequence.slice().reverse().map(base => complement[base]);
    return sequence.every((base, index) => base === reverseComplement[index]);
  }

  function findParameter(entries, key) {
    if (entries[key]) return [key, entries[key]];

    // RNA/RNA and modified RNA/RNA tables may be represented by the
    // reverse-complement-equivalent orientation.
    const [top, bottom] = key.split('/');
    const candidates = [
      `${reverseString(bottom)}/${reverseString(top)}`,
      `${bottom}/${top}`,
      `${reverseString(top)}/${reverseString(bottom)}`
    ];
    for (const candidate of candidates) {
      if (entries[candidate]) return [candidate, entries[candidate]];
    }
    return null;
  }

  function addRow(rows, label, parameter, values) {
    if (!values || !Number.isFinite(values.dH) || !Number.isFinite(values.dS) || !Number.isFinite(values.dG)) {
      throw new Error(`Invalid values for parameter “${parameter}”.`);
    }
    rows.push({ label, parameter, dH: values.dH, dS: values.dS, dG: values.dG });
  }

  function calculate() {
    try {
      if (!DB || !DB.sets) throw new Error('The parameter database was not loaded.');

      if ($('autofill').checked) generateComplement(false);

      const mode = $('mode').value;
      const condition = $('condition').value;
      const a = tokenize($('a').value);
      const b = tokenize($('b').value);

      if (a.length < 2) throw new Error('Enter a duplex containing at least two base pairs.');
      if (a.length !== b.length) throw new Error(`The two strands have different lengths (${a.length} and ${b.length} residues).`);
      if (!a.every(t => /^[ACGU]$/.test(t))) throw new Error('Strand 1 must be unmodified RNA (A/C/G/U).');
      if (!b.every(t => isExpectedStrand2Token(mode, t))) throw new Error('Strand 2 notation does not match the selected duplex type.');

      for (let i = 0; i < a.length; i++) {
        if (complement[canonicalBase(a[i])] !== canonicalBase(b[i])) {
          throw new Error(`The strands are not Watson–Crick complementary at position ${i + 1}.`);
        }
      }

      const setName = `${mode}_${condition}`;
      const set = DB.sets[setName];
      if (!set) throw new Error(`Parameter set not found: ${setName}.`);

      const rows = [];
      for (let i = 0; i < a.length - 1; i++) {
        const requestedKey = makeStepKey(mode, a.slice(i, i + 2), b.slice(i, i + 2));
        const found = findParameter(set.entries, requestedKey);
        if (!found) throw new Error(`No NN parameter was found for ${requestedKey} at step ${i + 1}.`);
        addRow(rows, `${i + 1}–${i + 2}`, found[0], found[1]);
      }

      const corrections = set.corrections || {};
      if (mode === 'RNA_DNA') {
        for (const position of [0, a.length - 1]) {
          const terminalBase = canonicalBase(a[position]);
          const correctionName = /[GC]/.test(terminalBase)
            ? 'init. rG−dC and rC−dG'
            : 'init. rA−dT or rU−dA';
          addRow(rows, position === 0 ? '5′ terminal initiation' : '3′ terminal initiation', correctionName, corrections[correctionName]);
        }
      } else {
        if (corrections.initiation) addRow(rows, 'initiation', 'initiation', corrections.initiation);
        for (const position of [0, a.length - 1]) {
          const terminalBase = canonicalBase(a[position]);
          const correctionName = mode === 'unmodified'
            ? (/[AU]/.test(terminalBase) ? 'per terminal AU' : 'per terminal GC')
            : (/[AU]/.test(terminalBase) ? 'per terminal AL or UH' : 'per terminal GJ or CK');
          addRow(rows, position === 0 ? '5′ terminal' : '3′ terminal', correctionName, corrections[correctionName]);
        }
      }

      const selfComplementary = isSelfComplementary(a);
      if (selfComplementary && corrections['Self-complementary']) {
        addRow(rows, 'symmetry', 'Self-complementary', corrections['Self-complementary']);
      }

      const totals = rows.reduce((sum, row) => ({
        dH: sum.dH + row.dH,
        dS: sum.dS + row.dS,
        dG: sum.dG + row.dG
      }), { dH: 0, dS: 0, dG: 0 });

      $('dh').textContent = totals.dH.toFixed(2);
      $('ds').textContent = totals.dS.toFixed(2);
      $('dg').textContent = totals.dG.toFixed(2);
      $('rows').innerHTML = rows.map(row =>
        `<tr><td>${escapeHtml(row.label)}</td><td>${escapeHtml(row.parameter)}</td><td>${row.dH.toFixed(2)}</td><td>${row.dS.toFixed(2)}</td><td>${row.dG.toFixed(2)}</td></tr>`
      ).join('');
      $('results').hidden = false;
      setMessage(`Calculation completed using ${setName}.${selfComplementary ? ' Self-complementary correction was applied automatically.' : ''}`, true);
      lastRows = rows;
    } catch (error) {
      $('results').hidden = true;
      setMessage(error.message, false);
    }
  }

  function escapeHtml(value) {
    return String(value).replace(/[&<>'"]/g, ch => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', "'": '&#39;', '"': '&quot;' }[ch]));
  }

  function downloadCsv() {
    if (!lastRows.length) {
      setMessage('Run a calculation before downloading CSV.', false);
      return;
    }
    const lines = [['Step', 'Parameter', 'dH', 'dS', 'dG'], ...lastRows.map(r => [r.label, r.parameter, r.dH, r.dS, r.dG])];
    const csv = lines.map(row => row.map(value => `"${String(value).replace(/"/g, '""')}"`).join(',')).join('\n');
    const link = document.createElement('a');
    link.href = URL.createObjectURL(new Blob([csv], { type: 'text/csv;charset=utf-8' }));
    link.download = 'NN_prediction.csv';
    link.click();
    URL.revokeObjectURL(link.href);
  }

  function notationText() {
    const mode = $('mode').value;
    const examples = {
      unmodified: 'Strand 2 notation: A, C, G, U',
      RNA_DNA: 'Strand 2 notation: dA, dC, dG, dT',
      OMe: 'Strand 2 notation: Ao, Co, Go, Uo',
      F: 'Strand 2 notation: Af, Cf, Gf, Uf',
      MOE: 'Strand 2 notation: Am, Cm, Gm, Um',
      LNA: 'Strand 2 notation: Al, Cl, Gl, Ul'
    };
    $('notation').textContent = examples[mode];
  }

  function setMessage(message, success) {
    $('msg').textContent = message;
    $('msg').className = success ? 'ok' : 'error';
  }

  function clearMessage() {
    $('msg').textContent = '';
    $('msg').className = '';
  }

  function onModeChange() {
    notationText();
    if ($('autofill').checked) generateComplement(false);
  }

  function init() {
    $('auto').addEventListener('click', () => generateComplement(true));
    $('calc').addEventListener('click', calculate);
    $('csv').addEventListener('click', downloadCsv);
    $('mode').addEventListener('change', onModeChange);
    $('a').addEventListener('input', () => {
      if ($('autofill').checked) generateComplement(false);
    });
    $('autofill').addEventListener('change', () => {
      if ($('autofill').checked) generateComplement(false);
    });
    notationText();
    if ($('autofill').checked) generateComplement(false);
    if (!DB) setMessage('The embedded parameter database could not be loaded.', false);
  }

  document.addEventListener('DOMContentLoaded', init);
})();
