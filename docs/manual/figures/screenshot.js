// Screenshot of an element of an HTML page, for the figures of the manual:
// by default the <pre> element (the pages of terminal.py).
//
// Usage: node screenshot.js PAGE.html OUT.png [SELECTOR [WIDTH]]
//   (requires Playwright); WIDTH: of the browser window (default 1280)
const { chromium } = require('playwright');
const path = require('path');

(async () => {
  const [page_path, out, selector = 'pre', width = '1280'] =
    process.argv.slice(2);
  const browser = await chromium.launch();
  const page = await browser.newPage({ deviceScaleFactor: 2,
    viewport: { width: parseInt(width), height: 800 } });
  await page.goto('file://' + path.resolve(page_path));
  await (await page.$(selector)).screenshot({ path: out });
  await browser.close();
})();
