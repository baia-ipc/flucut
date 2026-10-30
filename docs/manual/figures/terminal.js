// (c) 2026 Giorgio Gonnella
// A terminal transcript as an image, for the manual: the output of a
// program in a terminal (ANSI colors, e.g. the banner of flucut),
// captured with `script`, rendered in a terminal window. From the manual
// of Lumipore (github.com/baia-ipc/lumipore, docs/manual/figures), with
// the lines longer than the columns wrapped, as in a terminal.
//
// Usage: node terminal.js TRANSCRIPT OUTPUT.png COMMAND [COLUMNS]
//
// TRANSCRIPT: the output (ANSI escape sequences; carriage returns: only
// the text after the last one of each line is kept, as on the terminal);
// COMMAND: the command line shown after the prompt. Requires Playwright
// (with Chromium) and the fonts JetBrains Mono and FreeMono (Braille).
const { chromium } = require('playwright');
const fs = require('fs');

const BASIC = ['#1e1e1e', '#e5484d', '#2ebfa5', '#f5a524', '#5b9be6',
               '#c26dd1', '#3fb8c9', '#d4d4d4'];

function esc(s) {
  return s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
}

// the lines of the transcript as HTML (colored spans)
function render(text) {
  const out = [];
  let fg = null, bg = null, bold = false;
  for (let line of text.split('\n')) {
    line = line.replace(/\r+$/, '');
    const cr = line.lastIndexOf('\r');
    if (cr >= 0)
      line = line.slice(cr + 1);
    let html = '';
    for (const part of line.split(/(\x1b\[[0-9;]*m)/)) {
      const m = part.match(/^\x1b\[([0-9;]*)m$/);
      if (m) {
        const c = m[1] === '' ? [0] : m[1].split(';').map(Number);
        for (let i = 0; i < c.length; i++) {
          if (c[i] === 0) { fg = bg = null; bold = false; }
          else if (c[i] === 1) bold = true;
          else if (c[i] === 22) bold = false;
          else if (c[i] >= 30 && c[i] <= 37) fg = BASIC[c[i] - 30];
          else if (c[i] === 39) fg = null;
          else if (c[i] === 49) bg = null;
          else if ((c[i] === 38 || c[i] === 48) && c[i + 1] === 2) {
            const rgb = `rgb(${c[i + 2]},${c[i + 3]},${c[i + 4]})`;
            if (c[i] === 38) fg = rgb; else bg = rgb;
            i += 4;
          }
        }
      } else if (part) {
        const style = (fg ? `color:${fg};` : '') +
                      (bg ? `background:${bg};` : '') +
                      (bold ? 'font-weight:700;' : '');
        html += style ? `<span style="${style}">${esc(part)}</span>`
                      : esc(part);
      }
    }
    out.push(html);
  }
  while (out.length && out[out.length - 1] === '')
    out.pop();
  return out.join('\n');
}

(async () => {
  const [input, output, command, columns] = process.argv.slice(2);
  const cols = Number(columns || 100);
  const body = render(fs.readFileSync(input, 'utf8'));
  const page = `<!doctype html><html><head><meta charset="utf-8"><style>
body{margin:0;background:#fff;padding:12px}
.win{display:inline-block;background:#000;border-radius:10px;
box-shadow:0 4px 16px rgba(0,0,0,.25);overflow:hidden}
.bar{background:#2b2b2b;height:26px;display:flex;align-items:center;
padding:0 10px;gap:7px}
.dot{width:11px;height:11px;border-radius:50%}
pre{margin:0;padding:12px 14px 14px;color:#d4d4d4;
font:13px/1.18 'JetBrains Mono','FreeMono',monospace;width:${cols}ch;
font-variant-ligatures:none;font-feature-settings:"calt" 0,"liga" 0;
white-space:pre-wrap;word-break:break-all}
.prompt{color:#2ebfa5}
</style></head><body><div class="win"><div class="bar">
<span class="dot" style="background:#ff5f57"></span>
<span class="dot" style="background:#febc2e"></span>
<span class="dot" style="background:#28c840"></span></div>
<pre><span class="prompt">$</span> ${esc(command)}
${body}</pre></div></body></html>`;
  const browser = await chromium.launch();
  const p = await browser.newPage({ deviceScaleFactor: 2 });
  await p.setContent(page);
  await p.waitForTimeout(200);
  await p.locator('.win').screenshot({ path: output });
  await browser.close();
})();
