# RNA Duplex Stability Calculator

Open `index.html` directly in a browser or deploy the folder to GitHub Pages.
The parameter database is embedded in `js/parameters.js`, so the calculator works without a local web server.

## Input convention

- Strand 1: unmodified RNA, 5′→3′, using A/C/G/U.
- Strand 2: complementary strand, 3′→5′.
- Select **Automatically fill the Watson–Crick complementary strand** or click **Generate complement**.

Accepted Strand 2 notation:

- RNA: A/C/G/U
- DNA: dA/dC/dG/dT
- 2′-OMe: Ao/Co/Go/Uo
- 2′-F: Af/Cf/Gf/Uf
- 2′-MOE: Am/Cm/Gm/Um
- LNA: Al/Cl/Gl/Ul

## Updating parameters

Update `parameters.json`, then regenerate `js/parameters.js` as:

```javascript
window.NN_PARAMETER_DB = { ...contents of parameters.json... };
```


## Self-complementary correction

The symmetry correction is detected automatically. It is added only when the 5′→3′ sequence of Strand 1 is identical to its reverse complement, such as `AACCGGUU` for RNA (analogous to `AACCGGTT` for DNA). Merely entering a Watson–Crick complementary Strand 2 does not trigger this correction.
