// @ts-check
import { defineConfig } from 'astro/config';
import starlight from '@astrojs/starlight';
import starlightThemeNova from 'starlight-theme-nova';
import { readFileSync, existsSync } from 'node:fs';
import { fileURLToPath } from 'node:url';

const sidebarManifestPath = fileURLToPath(new URL('./build/sidebar.json', import.meta.url));
const referenceSidebar = existsSync(sidebarManifestPath)
	? JSON.parse(readFileSync(sidebarManifestPath, 'utf8'))
	: [];

// https://astro.build/config
export default defineConfig({
	integrations: [
		starlight({
			title: 'MIINT',
			description: 'Microbiome analytics for DuckDB.',
			plugins: [starlightThemeNova()],
			social: [
				{ icon: 'github', label: 'GitHub', href: 'https://github.com/the-miint/duckdb-miint' },
			],
			sidebar: [
				{
					label: 'Getting started',
					items: [
						{ label: 'Installation', slug: 'getting-started/installation' },
						{ label: 'Updating', slug: 'getting-started/updating' },
					],
				},
				{
					label: 'Internals',
					autogenerate: { directory: 'internals' },
				},
				{
					label: 'Reference',
					items: referenceSidebar,
				},
			],
		}),
	],
});
