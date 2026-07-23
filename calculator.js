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
      $('resultMeta').textContent = `${setName} · Strand A: ${$('a').value} · Strand B: ${$('b').value}${selfComplementary ? ' · self-complementary correction applied' : ''}`;
      setMessage(`Calculation completed using ${setName}.${selfComplementary ? ' Self-complementary correction was applied automatically.' : ''}`, true);
      lastRows = rows;
    } catch (error) {
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

  const translations = {
    en: {
      eyebrow:'WEB-BASED NEAREST-NEIGHBOR PREDICTION TOOL', title:'2′-Modified RNA Duplex Stability Calculator',
      intro:'Predict sequence-dependent Δ<i>H</i>°, Δ<i>S</i>°, and Δ<i>G</i>°<sub>37</sub> values for RNA/RNA, RNA/DNA, and four 2′-modified RNA/RNA duplexes under non-crowding and cell-like conditions.',
      duplexesLabel:'Duplexes:',duplexesText:'RNA/RNA · RNA/DNA · 2′-OMe RNA/RNA · 2′-F RNA/RNA · 2′-MOE RNA/RNA · LNA RNA/RNA',conditionsLabel:'Conditions:',conditionsText:'non-crowding · cell-like',outputsLabel:'Outputs:',outputsText:'Δ<i>H</i>° · Δ<i>S</i>° · Δ<i>G</i>°<sub>37</sub> · NN breakdown · CSV',language:'Language',howToUse:'How to use',step1Title:'Enter sequences',step1Text:'Input Strand A and auto-fill Strand B.',step2Title:'Select mode',step2Text:'Choose the duplex type and solution condition.',step3Title:'Calculate',step3Text:'Run the prediction and review the total values and NN breakdown.',sequenceInput:'Sequence Input',sequenceInputBold:'Sequence input.',sequenceInputHelp:'Enter Strand A and Strand B below. Strand B can be generated automatically.',acceptedNotationLabel:'Accepted notation:',acceptedNotation:'RNA, A/C/G/U; DNA, dA/dC/dG/dT; 2′-OMe, Ao/Co/Go/Uo; 2′-F, Af/Cf/Gf/Uf; 2′-MOE, Am/Cm/Gm/Um; LNA, Al/Cl/Gl/Ul.',strand1Label:'Strand A, 5′→3′ unmodified RNA strand',strand2Label:'Strand B, 3′→5′ complementary strand',autoFill:'Auto-fill Strand B',autoUpdate:'Update Strand B while Strand A or the mode changes',calculationMode:'Duplex type',condition:'Condition',modeRNA:'RNA/RNA',modeDNARNA:'RNA/DNA',modeOMe:'2′-OMe RNA/RNA',modeF:'2′-F RNA/RNA',modeMOE:'2′-MOE RNA/RNA',modeLNA:'LNA RNA/RNA',nonCrowding:'Non-crowding',cellLike:'Cell-like (40 wt% PEG 200)',modeSelectionBold:'Mode selection.',modeSelectionHelp:'Choose the duplex type and solution condition. The selected chemistry determines the required notation for Strand B.',selfNote:'The self-complementary correction is applied automatically only when a strand is identical to its own reverse-complement sequence.',calculate:'Calculate',downloadCsv:'Download CSV',executionBold:'Execution and output.',executionHelp:'Click Calculate to obtain Δ<i>H</i>°, Δ<i>S</i>°, Δ<i>G</i>°<sub>37</sub>, and the detailed NN contribution breakdown.',predictedParameters:'Predicted Thermodynamic Parameters',outputSubtitle:'Total values and NN contribution details',total:'total',nnContributions:'Nearest-neighbor contributions',stepCorrection:'Step/correction',parameter:'Parameter',parameterVersion:'Current parameter version:',references:'References',referenceNote:'The 2′-OMe, 2′-F, 2′-MOE, and LNA parameter sets used in the current website version are the parameter set dated 17 July 2026 and will be updated when the revised values become available.'
    },
    zh: {
      eyebrow:'基于网页的最近邻参数预测工具',title:'2′-修饰 RNA 双链稳定性计算器',intro:'预测 RNA/RNA、RNA/DNA 及四种 2′-修饰 RNA/RNA 双链在非拥挤和细胞样条件下的序列依赖性 Δ<i>H</i>°、Δ<i>S</i>° 和 Δ<i>G</i>°<sub>37</sub>。',duplexesLabel:'双链类型：',duplexesText:'RNA/RNA · RNA/DNA · 2′-OMe RNA/RNA · 2′-F RNA/RNA · 2′-MOE RNA/RNA · LNA RNA/RNA',conditionsLabel:'条件：',conditionsText:'非拥挤 · 细胞样',outputsLabel:'输出：',outputsText:'Δ<i>H</i>° · Δ<i>S</i>° · Δ<i>G</i>°<sub>37</sub> · NN 参数分解 · CSV',language:'语言',howToUse:'使用方法',step1Title:'输入序列',step1Text:'输入链 A，并自动生成链 B。',step2Title:'选择模式',step2Text:'选择双链类型和溶液条件。',step3Title:'计算',step3Text:'运行预测并查看总参数和 NN 分解结果。',sequenceInput:'序列输入',sequenceInputBold:'序列输入。',sequenceInputHelp:'在下方输入链 A 和链 B；链 B 可自动生成。',acceptedNotationLabel:'可用符号：',acceptedNotation:'RNA：A/C/G/U；DNA：dA/dC/dG/dT；2′-OMe：Ao/Co/Go/Uo；2′-F：Af/Cf/Gf/Uf；2′-MOE：Am/Cm/Gm/Um；LNA：Al/Cl/Gl/Ul。',strand1Label:'链 A：5′→3′ 未修饰 RNA 链',strand2Label:'链 B：3′→5′ 互补链',autoFill:'自动生成链 B',autoUpdate:'链 A 或模式变化时自动更新链 B',calculationMode:'双链类型',condition:'条件',modeRNA:'RNA/RNA',modeDNARNA:'RNA/DNA',modeOMe:'2′-OMe RNA/RNA',modeF:'2′-F RNA/RNA',modeMOE:'2′-MOE RNA/RNA',modeLNA:'LNA RNA/RNA',nonCrowding:'非拥挤',cellLike:'细胞样（40 wt% PEG 200）',modeSelectionBold:'模式选择。',modeSelectionHelp:'选择双链类型和溶液条件；所选化学类型决定链 B 的输入符号。',selfNote:'仅当单链与其反向互补序列完全相同时，程序才自动加入自互补修正。',calculate:'计算',downloadCsv:'下载 CSV',executionBold:'执行与输出。',executionHelp:'点击“计算”获得 Δ<i>H</i>°、Δ<i>S</i>°、Δ<i>G</i>°<sub>37</sub> 以及详细的 NN 参数贡献。',predictedParameters:'预测热力学参数',outputSubtitle:'总参数及 NN 贡献明细',total:'总值',nnContributions:'最近邻参数贡献',stepCorrection:'步骤/修正项',parameter:'参数',parameterVersion:'当前参数版本：',references:'参考文献',referenceNote:'当前网站使用的 2′-OMe、2′-F、2′-MOE 和 LNA 参数为 2026 年 7 月 17 日版本；修订参数确定后将进一步更新。'
    },
    ja: {
      eyebrow:'ウェブベース最近接塩基対予測ツール',title:'2′修飾 RNA 二重鎖安定性計算ツール',intro:'非クラウディングおよび細胞様条件における RNA/RNA、RNA/DNA、ならびに4種類の2′修飾 RNA/RNA 二重鎖の配列依存的な Δ<i>H</i>°、Δ<i>S</i>°、Δ<i>G</i>°<sub>37</sub> を予測します。',duplexesLabel:'二重鎖：',duplexesText:'RNA/RNA · RNA/DNA · 2′-OMe RNA/RNA · 2′-F RNA/RNA · 2′-MOE RNA/RNA · LNA RNA/RNA',conditionsLabel:'条件：',conditionsText:'非クラウディング · 細胞様',outputsLabel:'出力：',outputsText:'Δ<i>H</i>° · Δ<i>S</i>° · Δ<i>G</i>°<sub>37</sub> · NN 内訳 · CSV',language:'言語',howToUse:'使用方法',step1Title:'配列を入力',step1Text:'鎖 A を入力し、鎖 B を自動生成します。',step2Title:'モードを選択',step2Text:'二重鎖タイプと溶液条件を選択します。',step3Title:'計算',step3Text:'予測を実行し、合計値と NN 内訳を確認します。',sequenceInput:'配列入力',sequenceInputBold:'配列入力。',sequenceInputHelp:'鎖 A と鎖 B を入力してください。鎖 B は自動生成できます。',acceptedNotationLabel:'使用可能な表記：',acceptedNotation:'RNA：A/C/G/U；DNA：dA/dC/dG/dT；2′-OMe：Ao/Co/Go/Uo；2′-F：Af/Cf/Gf/Uf；2′-MOE：Am/Cm/Gm/Um；LNA：Al/Cl/Gl/Ul。',strand1Label:'鎖 A：5′→3′ 未修飾 RNA 鎖',strand2Label:'鎖 B：3′→5′ 相補鎖',autoFill:'鎖 B を自動入力',autoUpdate:'鎖 A またはモード変更時に鎖 B を更新',calculationMode:'二重鎖タイプ',condition:'条件',modeRNA:'RNA/RNA',modeDNARNA:'RNA/DNA',modeOMe:'2′-OMe RNA/RNA',modeF:'2′-F RNA/RNA',modeMOE:'2′-MOE RNA/RNA',modeLNA:'LNA RNA/RNA',nonCrowding:'非クラウディング',cellLike:'細胞様（40 wt% PEG 200）',modeSelectionBold:'モード選択。',modeSelectionHelp:'二重鎖タイプと溶液条件を選択してください。選択した化学種により鎖 B の表記が決まります。',selfNote:'自己相補性補正は、一本鎖が自身の逆相補配列と完全に一致する場合にのみ自動適用されます。',calculate:'計算',downloadCsv:'CSVをダウンロード',executionBold:'計算と出力。',executionHelp:'「計算」をクリックすると、Δ<i>H</i>°、Δ<i>S</i>°、Δ<i>G</i>°<sub>37</sub> と NN 寄与の詳細が表示されます。',predictedParameters:'予測熱力学パラメータ',outputSubtitle:'合計値および NN 寄与の詳細',total:'合計',nnContributions:'最近接塩基対パラメータの寄与',stepCorrection:'ステップ/補正',parameter:'パラメータ',parameterVersion:'現在のパラメータ版：',references:'参考文献',referenceNote:'現在のウェブサイトで使用している 2′-OMe、2′-F、2′-MOE、および LNA パラメータは2026年7月17日版であり、改訂値が確定後に更新されます。'
    },
    ko: {
      eyebrow:'웹 기반 최근접 이웃 예측 도구',title:'2′-변형 RNA 이중가닥 안정성 계산기',intro:'비혼잡 및 세포 유사 조건에서 RNA/RNA, RNA/DNA 및 네 종류의 2′-변형 RNA/RNA 이중가닥에 대한 서열 의존적 Δ<i>H</i>°, Δ<i>S</i>°, Δ<i>G</i>°<sub>37</sub> 값을 예측합니다.',duplexesLabel:'이중가닥:',duplexesText:'RNA/RNA · RNA/DNA · 2′-OMe RNA/RNA · 2′-F RNA/RNA · 2′-MOE RNA/RNA · LNA RNA/RNA',conditionsLabel:'조건:',conditionsText:'비혼잡 · 세포 유사',outputsLabel:'출력:',outputsText:'Δ<i>H</i>° · Δ<i>S</i>° · Δ<i>G</i>°<sub>37</sub> · NN 세부내역 · CSV',language:'언어',howToUse:'사용 방법',step1Title:'서열 입력',step1Text:'가닥 A를 입력하고 가닥 B를 자동 생성합니다.',step2Title:'모드 선택',step2Text:'이중가닥 유형과 용액 조건을 선택합니다.',step3Title:'계산',step3Text:'예측을 실행하고 총값과 NN 세부내역을 확인합니다.',sequenceInput:'서열 입력',sequenceInputBold:'서열 입력.',sequenceInputHelp:'아래에 가닥 A와 가닥 B를 입력하십시오. 가닥 B는 자동 생성할 수 있습니다.',acceptedNotationLabel:'허용 표기:',acceptedNotation:'RNA: A/C/G/U; DNA: dA/dC/dG/dT; 2′-OMe: Ao/Co/Go/Uo; 2′-F: Af/Cf/Gf/Uf; 2′-MOE: Am/Cm/Gm/Um; LNA: Al/Cl/Gl/Ul.',strand1Label:'가닥 A: 5′→3′ 비변형 RNA 가닥',strand2Label:'가닥 B: 3′→5′ 상보 가닥',autoFill:'가닥 B 자동 입력',autoUpdate:'가닥 A 또는 모드 변경 시 가닥 B 업데이트',calculationMode:'이중가닥 유형',condition:'조건',modeRNA:'RNA/RNA',modeDNARNA:'RNA/DNA',modeOMe:'2′-OMe RNA/RNA',modeF:'2′-F RNA/RNA',modeMOE:'2′-MOE RNA/RNA',modeLNA:'LNA RNA/RNA',nonCrowding:'비혼잡',cellLike:'세포 유사 (40 wt% PEG 200)',modeSelectionBold:'모드 선택.',modeSelectionHelp:'이중가닥 유형과 용액 조건을 선택하십시오. 선택한 화학 유형에 따라 가닥 B의 표기가 결정됩니다.',selfNote:'자기상보성 보정은 한 가닥이 자신의 역상보 서열과 완전히 동일한 경우에만 자동 적용됩니다.',calculate:'계산',downloadCsv:'CSV 다운로드',executionBold:'계산 및 출력.',executionHelp:'계산을 클릭하면 Δ<i>H</i>°, Δ<i>S</i>°, Δ<i>G</i>°<sub>37</sub> 및 상세 NN 기여도를 확인할 수 있습니다.',predictedParameters:'예측 열역학 파라미터',outputSubtitle:'총값 및 NN 기여도 상세',total:'총값',nnContributions:'최근접 이웃 기여도',stepCorrection:'단계/보정',parameter:'파라미터',parameterVersion:'현재 파라미터 버전:',references:'참고문헌',referenceNote:'현재 웹사이트에 사용된 2′-OMe, 2′-F, 2′-MOE 및 LNA 파라미터는 2026년 7월 17일 버전이며, 개정 값이 확정되면 업데이트됩니다.'
    }
  };

  let currentLanguage = 'en';
  function t(key){ return (translations[currentLanguage] && translations[currentLanguage][key]) || translations.en[key] || key; }
  function applyLanguage(lang){
    currentLanguage = translations[lang] ? lang : 'en';
    document.documentElement.lang = currentLanguage === 'zh' ? 'zh-CN' : currentLanguage;
    document.querySelectorAll('[data-i18n]').forEach(el => { const value=t(el.dataset.i18n); if(value != null) el.textContent=value; });
    document.querySelectorAll('[data-i18n-html]').forEach(el => { const value=t(el.dataset.i18nHtml); if(value != null) el.innerHTML=value; });
    notationText();
    if (lastRows.length) calculate();
  }

  function notationText() {
    const mode = $('mode').value;
    const examples = {
      en:{unmodified:'Strand B notation: A, C, G, U',RNA_DNA:'Strand B notation: dA, dC, dG, dT',OMe:'Strand B notation: Ao, Co, Go, Uo',F:'Strand B notation: Af, Cf, Gf, Uf',MOE:'Strand B notation: Am, Cm, Gm, Um',LNA:'Strand B notation: Al, Cl, Gl, Ul'},
      zh:{unmodified:'链 B 符号：A、C、G、U',RNA_DNA:'链 B 符号：dA、dC、dG、dT',OMe:'链 B 符号：Ao、Co、Go、Uo',F:'链 B 符号：Af、Cf、Gf、Uf',MOE:'链 B 符号：Am、Cm、Gm、Um',LNA:'链 B 符号：Al、Cl、Gl、Ul'},
      ja:{unmodified:'鎖 B の表記：A、C、G、U',RNA_DNA:'鎖 B の表記：dA、dC、dG、dT',OMe:'鎖 B の表記：Ao、Co、Go、Uo',F:'鎖 B の表記：Af、Cf、Gf、Uf',MOE:'鎖 B の表記：Am、Cm、Gm、Um',LNA:'鎖 B の表記：Al、Cl、Gl、Ul'},
      ko:{unmodified:'가닥 B 표기: A, C, G, U',RNA_DNA:'가닥 B 표기: dA, dC, dG, dT',OMe:'가닥 B 표기: Ao, Co, Go, Uo',F:'가닥 B 표기: Af, Cf, Gf, Uf',MOE:'가닥 B 표기: Am, Cm, Gm, Um',LNA:'가닥 B 표기: Al, Cl, Gl, Ul'}
    };
    $('notation').textContent = examples[currentLanguage][mode];
  }

  function setMessage(message, success) { $('msg').textContent = message; $('msg').className = success ? 'ok' : 'error'; }
  function clearMessage() { $('msg').textContent = ''; $('msg').className = ''; }
  function onModeChange() { notationText(); if ($('autofill').checked) generateComplement(false); calculate(); }

  function init() {
    $('auto').addEventListener('click', () => { generateComplement(true); calculate(); });
    $('calc').addEventListener('click', calculate);
    $('csv').addEventListener('click', downloadCsv);
    $('mode').addEventListener('change', onModeChange);
    $('condition').addEventListener('change', calculate);
    $('language').addEventListener('change', e => applyLanguage(e.target.value));
    $('a').addEventListener('input', () => { if ($('autofill').checked) generateComplement(false); });
    $('autofill').addEventListener('change', () => { if ($('autofill').checked) generateComplement(false); });
    applyLanguage('en');
    generateComplement(false);
    calculate();
  }

  document.addEventListener('DOMContentLoaded', init);
})();
