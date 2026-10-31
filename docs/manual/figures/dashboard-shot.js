// A full-page screenshot of the live-mode dashboard, at a wide viewport.
//
// Usage: node dashboard-shot.js <flucut_report.html> <out.png>
//
// Requires Playwright with Chromium (as screenshot.js).
const { chromium } = require('playwright');
const path = require('path');

(async () => {
  const [html, out] = process.argv.slice(2);
  const browser = await chromium.launch();
  const page = await browser.newPage({
    viewport: { width: 1680, height: 945 }, deviceScaleFactor: 2 });
  await page.goto('file://' + path.resolve(html));
  await page.waitForTimeout(600);
  await page.screenshot({ path: out });
  await browser.close();
})();
