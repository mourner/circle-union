import config from 'eslint-config-mourner';

// Skip generated artifacts: the UMD bundle (`npm run build`) and emitted type declarations.
export default [{ignores: ['circle-union.js', '**/*.d.ts']}, ...config];
